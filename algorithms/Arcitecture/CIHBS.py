import math

from Bio.PDB import Selection, Atom, NeighborSearch

from Arcitecture import BaseEnums

class CIHBS:
    def __init__(self):
        self.newCIHBS = []
        self.mapDonorAcseptor = []

    def setParameters(self, structure):
        self.alphaCarbonChain = [atom for atom in Selection.unfold_entities(structure, 'A') if atom.get_name() == "CA"]
        self.maimChain = [atom for atom in Selection.unfold_entities(structure, 'A') if
                          (atom.get_name() == "CA" or len(atom.get_name()) == 1)]
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    def getCarbonChain(self):
        return self.alphaCarbonChain

    def setNeighbourSearch(self, structure):
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    def findTargetNeigbourExept(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if(elem.element in target and elem != atom and elem.get_parent() != atom.get_parent())]
        return neighbours

    def findTargetNeigbour(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if elem.element in target]
        return neighbours

    def findNeighboursExeptAtom(self, atom, radius, condotion):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if (elem.element != condotion or elem == atom)]
        return neighbours

    def connectAtoms(self, atom1, atom2, bond):
        #print("connecting", atom1, atom2)
        self.newCIHBS.append([atom1, atom2, bond])

    def checkInnerGroups(self, structure):

        cihbsList = []
        for residue in structure.get_residues():
            #print("\nfor residue = ", residue.get_resname())

            # проверка на бензольные кольца и их добавление
            if (residue.get_resname() in BaseEnums.Groups.sycleResiduces.value):
                cycle = list(residue.get_atoms())
                atomsName = [x.get_name() for x in cycle]

                try:
                    if(residue.get_resname() == "PRO"):
                        #print("adding CA-N for PRO", cycle[:atomsName.index("C")])
                        cihbsList.append(cycle[:atomsName.index("C")])
                        cycle = cycle[atomsName.index("CB"):]
                    else:
                        cycle = cycle[atomsName.index("CG"):]
                except ValueError:
                    #print("Не найден опорный атом CG")
                    pass

                if(residue.get_resname() == "PRO"):
                    cycle = [elem for elem in cycle if(elem.element == "C" and elem.get_name() != "CA")]

                #print("atoms cycle", cycle)
                del atomsName
                cihbsList.append(cycle)

            for atom in residue:
                #print("for atom ", atom.get_name(), type(atom.get_name()))

                # проверка простых групп
                if (residue.get_resname() in BaseEnums.Groups.singleResidues.value):
                    if (len(atom.get_name()) > 1 and atom.get_name()[0] in BaseEnums.Atoms.singleTargetElements.value):
                        self.setNeighbourSearch(residue)
                        neighbour = self.findTargetNeigbour(atom, 2.0, BaseEnums.Atoms.allAtomsElements.value)
                        #print("Neighbours = ", len(neighbour), neighbour)

                        if (len(neighbour) > 2):
                            neighbour.pop()
                            #print("Превышение количества соседей для ССИВС в простой группе")

                        #print("Neighbours simple group = ", len(neighbour), neighbour)
                        cihbsList.append(neighbour)

                # проверка резонансных групп, проверка связей N-C-O, O-C-O, N-C-N и тупиковые C-O
                if (atom.element == "C"):
                    self.setNeighbourSearch(structure)
                    neighbour = self.findNeighboursExeptAtom(atom, 2.0, "C")
                    atomO = Atom.Atom("O", None, None, None, '', "O", None, "O")

                    if (len(neighbour) > 3 or atomO.element in [elem.element for elem in neighbour]):
                        #print("Neighbours resonance = ", len(neighbour), neighbour)
                        cihbsList.append(neighbour)

        self.innerCIHBS = cihbsList
        return cihbsList

    def removeCarbon(self, group):
        tmp = []
        for elem in group:
            #print("elem.element in removeCarbon", elem.element)
            if(elem.get_name()[0] != "C"):
                tmp.append(elem)
        return tmp

    def getResidueAtoms(self, atom):
        #print("parent = ", atom.get_parent())
        return list(atom.get_parent().get_atoms())

    def getResseq(self, residue):
        resseq = [st.split('=') for st in str(residue).split(' ')]
        return resseq[4][1]

    def findInnerCIHBS(self, residue):

        for group in self.innerCIHBS:
            tmp = [self.getResseq(i.get_parent()) for i in group]

            if(len(tmp) == tmp.count(str(self.getResseq(residue))) and tmp[0] == str(self.getResseq(residue)) ):
                #print("TMP = ", tmp, type(tmp), str(self.getResseq(residue)), tmp[0] == str(self.getResseq(residue)) )
                group = self.removeCarbon(list(group))
                return group

        #print("is it GLU? LEU? VAL?")
        return None

    def findClosest(self, atoms, target):
        distances = [math.dist(i.get_coord(), target.get_coord()) for i in atoms]
        #print("findClosest to", atoms[distances.index(min(distances))],
         #     atoms[distances.index(min(distances))].get_parent())
        return atoms[distances.index(min(distances))]

    def findAlternativePhysicalOperator(self, target_element, find_element):
        self.setNeighbourSearch(sum(self.innerCIHBS, []))
        neighbours = [i for i in self.findTargetNeigbour(target_element, 4.0, find_element)
                      if len(i.get_name()) > 1]
        #print("chain connection failed", str(target_element), target_element.get_parent())
        if (neighbours != []):
            #print("no N found, found alternative variant",
                  #[i.get_parent() for i in neighbours], neighbours)
            #print("closest is", self.findClosest(neighbours, target_element).get_parent(),
                  #self.findClosest(neighbours, target_element))
            self.connectAtoms(self.findClosest(neighbours, target_element), target_element,
                              BaseEnums.CHIBSBond.physicalOperator)

    def connectPhysicalOperators(self):

        # получаем пентофрагмент
        for i in range(5, len(self.alphaCarbonChain)):
            group = self.alphaCarbonChain[i-5:i]

            self.setNeighbourSearch(self.getResidueAtoms(group[0]))
            i_4_element = self.findTargetNeigbour(group[0], 4.0, "O")[0]
            #i_3_element = self.findTargetNeigbour(self.getResidueAtoms(group[1]), group[1], 4.0, "O")[0]
            self.setNeighbourSearch(self.getResidueAtoms(group[4]))
            i_0_element = self.findTargetNeigbour(group[4], 2.0, "N")[0]
            #print(i_4_element, i_0_element)

            if(group[4] not in BaseEnums.Groups.antiConnectors.value):
                # получаем внутренние ССИВС для i-го элемента
                #print("for alpha group", [j.get_parent() for j in group])
                groupCIHBS = self.findInnerCIHBS(group[4].get_parent())
                if(groupCIHBS):
                    #print("groupCIHBS", [i.get_parent() for i in groupCIHBS], groupCIHBS)
                    pass

                #если SER, THR перекрывают, то соединяем остаток и 3
                #if(groupCIHBS and len(groupCIHBS) == 1 and math.dist(groupCIHBS[0].get_coord(), i_0_element.get_coord()) <= 4):
                #    print("connecting in res SER THR = ", math.dist(groupCIHBS[0].get_coord(), i_0_element.get_vector()), groupCIHBS[0].get_parent() )
                #    self.connectAtoms(i_3_element, groupCIHBS[0], BaseEnums.CHIBSBond.outerCIHBS)

                #иначе соединяем 0 и 4
                if(math.dist(i_4_element.get_coord(), i_0_element.get_coord()) <= 4):
                    #print("connecting in chain = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))
                    self.connectAtoms(i_0_element, i_4_element, BaseEnums.CHIBSBond.physicalOperator)
                else:
                    self.findAlternativePhysicalOperator(i_4_element, 'N')
                    self.findAlternativePhysicalOperator(i_0_element, 'O')

                # соединяем с остатком
                # сделать проверку на двойную связь!
                # + проверка на пересечение SER THR
                #if groupCIHBS:

                if(groupCIHBS and math.dist(self.findClosest(groupCIHBS, i_4_element).get_vector(), i_4_element.get_vector()) <= 4):
                    #print("connecting residual_4", [i.get_parent() for i in groupCIHBS], "with DIST = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))
                    self.connectAtoms(self.findClosest(groupCIHBS, i_4_element), i_4_element, BaseEnums.CHIBSBond.physicalOperator)

                #print("\n")
        #print(self.newCIHBS)
        return self.newCIHBS

    def removeNitrogen(self, group):
        tmp = []
        for elem in group:
            #print("elem.element in removeCarbon", elem.element)
            if( (elem.get_name()[0] == "N" and len(elem.get_name())) != 1 or elem.get_name()[0] != "N"):
                tmp.append(elem)
        return tmp

    def cleanMainChainAndCarbon(self, group):
        tmp = []
        for elem in group:
            if ( len(elem.get_name()) != 1 and elem.get_name()[0] != "C"):
                tmp.append(elem)
        return tmp

    def checkForIons(self, group, Ion):
        if(Ion == BaseEnums.Atoms.donor):
            group = [i for i in group if( i.get_name() not in BaseEnums.Groups.plusIONResidue.value[1]
                                           and i.get_parent() not in BaseEnums.Groups.plusIONResidue.value[0] )]
        else:
            group = [i for i in group if ( i.get_name() not in BaseEnums.Groups.minusIONResidue.value[1]
                                           and i.get_parent() not in BaseEnums.Groups.minusIONResidue.value[0] )]
        #print("checkForIons after cut", group)
        return group

    def checkForMainChain(self, group):
        return [i for i in group if( i not in self.maimChain)]

    def sordByDistance(self, target, group):
        pass

    def connectResidueCIHBS(self):

        #чистим
        cleanedCIHBS = []
        map = []

        for group in self.innerCIHBS:
            cleanedCIHBS.append(self.cleanMainChainAndCarbon(group))
            map.append(self.removeNitrogen(self.removeCarbon(group)))
        #print("cleanedCIHBS", cleanedCIHBS)

        gotInitiatorAcid = False
        for group in cleanedCIHBS:
            #print("for res", [i.get_parent() for i in group], group)

            for elem in group:
                #if(elem.get_parent().get_resname() in BaseEnums.Groups.CIHBSInitialisators.value):
                #    gotInitiatorAcid = True

                self.setNeighbourSearch( sum(map, []) )
                neighbours = self.findTargetNeigbourExept(elem, 3.7, ["N", "O", "S"])
                #print("neighbours", [i.get_parent() for i in neighbours], neighbours)


                #проверяем на ионы
                neighbours = self.checkForIons(neighbours, BaseEnums.Atoms.donor
                                    if(elem.element == BaseEnums.Atoms.donor.value) else BaseEnums.Atoms.acceptor)
                #проверяем на входжение в основную цепь
                #neighbours = self.checkForMainChain(neighbours)
                """if (neighbours != []):
                    neighbours = self.findClosest(neighbours, elem)
                    print("findClosest: ", neighbours)"""
                if (neighbours != []):
                    for neighbour in neighbours:
                        self.connectAtoms(neighbour, elem, BaseEnums.CHIBSBond.acceptorAcceptor if(elem.element == neighbour.element and elem.element == "O") else BaseEnums.CHIBSBond.acceptorAcceptor)



            #print("\n")
        return self.newCIHBS




