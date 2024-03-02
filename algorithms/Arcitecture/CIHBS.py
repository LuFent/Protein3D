
import math

from Bio.PDB import Selection, Atom, NeighborSearch
import warnings

from . import BaseEnums


class CIHBS:
    def __init__(self):
        self.newCIHBS = []

    def setParameters(self, structure):
        self.alphaCarbonChain = [atom for atom in Selection.unfold_entities(structure, 'A') if atom.get_name() == "CA"]
        self.mainChain = [atom for atom in Selection.unfold_entities(structure, 'A') if
                          ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                           and atom.get_parent().get_resname() != "HOH")]
        self.mainChainDict = {atom.get_parent(): atom.get_name() for atom in Selection.unfold_entities(structure, 'A') if
                          ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                           and atom.get_parent().get_resname() != "HOH")}
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

        self.nsStructure = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    def getMainChain(self):
        return self.mainChain

    def getCarbonChain(self):
        return self.alphaCarbonChain

    def getInnerCIHBS(self):
        if self.innerCIHBS:
            return self.innerCIHBS
        raise ValueError("Возвращается пустой объект.")

    def getNewCIHBS(self):
        if self.newCIHBS:
            return self.newCIHBS
        raise ValueError("Возвращается пустой объект.")

    def setNeighbourSearch(self, structure):
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    def findTargetNeigbourExept(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if
                      (elem.element in target and elem != atom and elem.get_parent() != atom.get_parent())]
        return neighbours

    def findTargetNeigbour(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if elem.element in target]
        return neighbours

    def connectAtoms(self, atom1, atom2, bond):
        # print("connecting", atom1, atom2)
        self.newCIHBS.append([atom1, atom2, bond])



    def checkInnerGroups(self, structure):

        cihbsList = []
        for residue in structure.get_residues():
            if residue.get_resname() != "HOH":

                #print("\nfor residue = ", residue)
                # checking for benzene rings and adding them
                if residue.get_resname() in BaseEnums.Groups.sycleResidues.value:
                    cycle = list(residue.get_atoms())
                    atomsName = [x.get_name() for x in cycle]

                    try:
                        cycle = cycle[atomsName.index("CG"):]
                        if residue.get_resname() == "PRO":
                            cycle = [elem for elem in cycle if (elem.get_name() not in ["C", "O"])]
                    except ValueError:
                        warnings.warn("Не найден опорный атом CG для {}. Проверьте целостность данных.".format(residue))

                    # print("atoms cycle", cycle)
                    del atomsName
                    cihbsList.append(cycle)

                try:
                    self.setNeighbourSearch(residue)
                    # simple residues CN, CS, CO and etc...
                    if residue.get_resname() in BaseEnums.Groups.singleResidues.value.keys():
                        try:
                            atom = BaseEnums.Groups.singleResidues.value.get(residue.get_resname())
                            atom_obj = list(residue.get_atoms())[
                                [element.get_name() for element in residue.get_atoms()].index(atom)]
                            neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)

                            cihbsList.append(neighbour)
                        except ValueError:
                            warnings.warn("Не найден опорный атом в двуплетах {}".format(residue))

                    # all troplets OCO, NCN and etc...
                    if residue.get_resname() in BaseEnums.Groups.complexResidues.value.keys():
                        try:
                            atom = BaseEnums.Groups.complexResidues.value.get(residue.get_resname())
                            atom_obj = list(residue.get_atoms())[
                                [element.get_name() for element in residue.get_atoms()].index(atom)]
                            neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)
                            neighbour = neighbour[[element.get_name() for element in neighbour].index(atom):]

                            # print("Neighbours complex group = ", len(neighbour), neighbour)
                            cihbsList.append(neighbour)
                        except ValueError:
                            warnings.warn("Не найден опорный атом в триплетах {}".format(residue) )
                            pass
                except IndexError:
                    warnings.warn("Ошибка распаковки аминокислоты {}. Проверьте целостность данных.".format(residue))


        # NCO
        tmp = [i for i in self.mainChain if (i.get_name() != "CA")]
        cihbsList.append(tmp)

        cihbsList = sum(cihbsList, [])
        cihbsDict = {}
        for i in cihbsList:
            if i.get_parent() in cihbsDict.keys():
                tmp_list = cihbsDict.get(i.get_parent())
                tmp_list.append(i)
                cihbsDict.update({i.get_parent(): tmp_list})
                del tmp_list
            else:
                cihbsDict[i.get_parent()] = [i]

        self.innerCIHBSDict = cihbsDict
        self.innerCIHBS = cihbsList



    def getResseq(self, residue):
        resseq = [st.split('=') for st in str(residue).split(' ')]
        return resseq[4][1]

    def findClosest(self, atoms, target):
        atoms = [i for i in atoms if (i.element != "C"
                                      and i.get_parent() not in BaseEnums.Groups.doubleBondAcid.value.keys()
                                      and i.get_name() not in BaseEnums.Groups.doubleBondAcid.value.values()
                                      and not i.is_disordered())]
        try:
            distances = [math.dist(i.get_coord(), target.get_coord()) for i in atoms]
            return atoms[distances.index(min(distances))]
        except ValueError:
            return None

    def findAlternativePhysicalOperator(self, target_element, find_element):

        neighbours = [i for i in self.findTargetNeigbour(target_element, 4.0, find_element)]
        # if abs( int(self.getResseq(i.get_parent())) - int(self.getResseq(target_element.get_parent())) ) > 4]
        print("chain connection failed", str(target_element), target_element.get_parent())
        if neighbours != [] and neighbours:
            print("no ATOM found, found alternative variant",
                  [i.get_parent() for i in neighbours], neighbours)
            try:
                print("closest is", self.findClosest(neighbours, target_element).get_parent(),
                      self.findClosest(neighbours, target_element))

                self.connectAtoms(self.findClosest(neighbours, target_element), target_element,
                                  BaseEnums.CHIBSBond.physicalOperator)
            except AttributeError:
                print("nobody nearby...")

    def connectPhysicalOperators(self):
        self.setNeighbourSearch(self.innerCIHBS)

        # получаем пентофрагмент
        for i in range(5, len(self.alphaCarbonChain)):
            group_4 = self.alphaCarbonChain[i - 5].get_parent()
            group_0 = self.alphaCarbonChain[i].get_parent()

            # if group[4].get_parent().get_resname() not in BaseEnums.Groups.antiConnectors.value:
            # получаем внутренние ССИВС для i-го элемента
            #print("for alpha group", [j.get_parent() for j in group])

            groupCIHBS = self.innerCIHBSDict.get(group_4)
            #print("groupCIHBS", [i.get_parent() for i in groupCIHBS], groupCIHBS)

            def getTargetElementInResidue(group, target):
                return  list(group.get_atoms())[list(i.get_name() for i in group.get_atoms()).index(target)]
            i_0_element = getTargetElementInResidue(group_4, "N")
            i_4_element = getTargetElementInResidue(group_0, "O")

            groupCIHBS = groupCIHBS[:groupCIHBS.index(i_0_element)]

            # если SER, THR перекрывают, то соединяем остаток и 3
            # if(groupCIHBS and len(groupCIHBS) == 1 and math.dist(groupCIHBS[0].get_coord(), i_0_element.get_coord()) <= 4):
            #    print("connecting in res SER THR = ", math.dist(groupCIHBS[0].get_coord(), i_0_element.get_vector()), groupCIHBS[0].get_parent() )
            #    self.connectAtoms(i_3_element, groupCIHBS[0], BaseEnums.CHIBSBond.outerCIHBS)

            # иначе соединяем 0 и 4
            if math.dist(i_4_element.get_coord(), i_0_element.get_coord()) <= 4:
                #print("connecting in chain = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))
                self.connectAtoms(i_0_element, i_4_element, BaseEnums.CHIBSBond.physicalOperator)
            else:
                pass
                # self.findAlternativePhysicalOperator(i_4_element, 'N')
                # !!!!self.findAlternativePhysicalOperator(i_0_element, 'O')

            # соединяем с остатком
            # сделать проверку на двойную связь
            closest = self.findClosest(groupCIHBS, i_4_element)

            if groupCIHBS and closest and math.dist(closest.get_vector(), i_4_element.get_vector()) <= 4 \
                    and group_4.get_resname() not in BaseEnums.Groups.antiConnectors.value:
                #print("connecting residual_4 to ", closest, "with DIST = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))

                self.connectAtoms(closest, i_4_element, BaseEnums.CHIBSBond.acceptorAcceptor
                if (i_4_element.element == closest.element)
                else BaseEnums.CHIBSBond.acceptorAcceptor)

            #print("\n")


    def checkForIons(self, group, Ion):
        if (Ion == BaseEnums.Atoms.donor):
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.plusIONResidue.value.keys()
                                          and i.get_parent() not in BaseEnums.Groups.plusIONResidue.value.values()]
        else:
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.minusIONResidue.value.keys()
                                          and i.get_parent() not in BaseEnums.Groups.minusIONResidue.value.values()]
        # print("checkForIons after cut", group)
        return group

    def checkRetrospective(self, residue_resseq, neighbours):
        return [i for i in neighbours if int(self.getResseq(i.get_parent())) < int(residue_resseq)]

    def connectResidueCIHBS(self):

        # чистим
        cleanedCIHBS = [i for i in self.innerCIHBS
                        if i.element != "C"
                        and i.get_name() != "N"]

        self.setNeighbourSearch(cleanedCIHBS)
        for elem in cleanedCIHBS:
            print("for element", elem, elem.get_parent())

            neighbours = self.findTargetNeigbourExept(elem, 3.7, ["N", "O", "S"])
            print("neighbours", [i.get_parent() for i in neighbours], neighbours)

            # проверяем на ионы
            neighbours = self.checkForIons(neighbours, BaseEnums.Atoms.donor
            if (elem.element == BaseEnums.Atoms.donor.value) else BaseEnums.Atoms.acceptor)

            # проверяем на ретроспективность ССИВС
            neighbours = self.checkRetrospective(self.getResseq(elem.get_parent()), neighbours)
            print("neighbours after sort", [i.get_parent() for i in neighbours], neighbours)

            if (neighbours != []):
                for neighbour in neighbours:
                    self.connectAtoms(neighbour, elem,
                                      BaseEnums.CHIBSBond.acceptorAcceptor if (
                                                  elem.element == neighbour.element and elem.element == "O")
                                      else BaseEnums.CHIBSBond.acceptorAcceptor)

            print("\n")

