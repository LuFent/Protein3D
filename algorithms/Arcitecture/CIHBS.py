import math

from Bio.PDB import Selection, NeighborSearch

from . import BaseEnums


class CIHBS:
    def __init__(self):
        self.newCIHBS = []

    def setParameters(self, structure):
        self.alphaCarbonChain = [atom for atom in Selection.unfold_entities(structure, 'A')
                                 if atom.get_name() == "CA"]
        self.mainChain = [atom for atom in Selection.unfold_entities(structure, 'A')
                          if ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                              and atom.get_parent().get_resname() != "HOH")]
        self.mainChainDict = {atom.get_parent(): atom.get_name() for atom in Selection.unfold_entities(structure, 'A')
                              if ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                                  and atom.get_parent().get_resname() != "HOH")}

        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))
        self.nsStructure = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    # получить атомы основной цепи
    def getMainChain(self):
        return [i.get_serial_number() for i in self.mainChain]

    # получить альфауглерод
    def getCarbonChain(self):
        return self.alphaCarbonChain

    # получить номера атомов простых ССИВС
    def getInnerCIHBS(self):
        if self.innerCIHBS:
            return [atom.get_serial_number() for atom in self.innerCIHBS]
        raise ValueError("Возвращается пустой объект.")

    # получить сложные ССИВС между аминокислотами
    def getNewCIHBS(self):
        if self.newCIHBS:
            dictinary = {'Hydrogen': [self.getR_O(), self.getR_R(), self.getOthers()],
                         'Physics': self.getNH_O()}
            return dictinary
        raise ValueError("Возвращается пустой объект.")

    # получение физических операторов между i и i-n элементами пентофрагментов
    def getNH_O(self):
        return {'NH_Oi_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 3
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 4
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 5
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 6
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O']}

    # получение ССИВС между iм остатком и i-n кислородом пентофрагментов
    def getR_O(self):
        return {'Ri_Oi_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 3
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 4
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 5
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 6
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O']}

    # получение ССИВС между iм остатком и i-n остатком
    def getR_R(self):
        return {'Ri_Ri_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 3
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 4
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 5
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                            if
                            abs(int(self.getResseq(i[0].get_parent())) - int(self.getResseq(i[1].get_parent()))) == 6
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2]}

    # получение остальных ССИВС
    def getOthers(self):
        return {'Others': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.newCIHBS
                           if
                           i[2] != BaseEnums.CHIBSBond.physicalOperator
                           and ((len(i[0].get_name()) >= 2 and i[1].get_name() != 'O')
                                or (len(i[0].get_name()) == 1 and len(i[1].get_name()) >= 2))]}

    # устанавливаем маску поиска по атомам
    def setNeighbourSearch(self, structure):
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    # поиск соседних целевых атомов среди соседних аминокислот
    def findTargetNeigbourExept(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if
                      (elem.element in target and elem != atom and elem.get_parent() != atom.get_parent())]
        return neighbours

    # поиск соседних целевых атомов
    def findTargetNeigbour(self, atom, radius, target):
        neighbours = [elem for elem in self.ns.search(atom.get_coord(), radius, 'A') if elem.element in target]
        return neighbours

    # соединяем атомы
    def connectAtoms(self, atom1, atom2, bond):
        self.newCIHBS.append([atom1, atom2, bond])

    def checkInnerGroups(self, structure):

        cihbsList = []
        for residue in structure.get_residues():
            if residue.get_resname() != "HOH":

                # циклические аминокислоты
                if residue.get_resname() in BaseEnums.Groups.sycleResidues.value:
                    cycle = list(residue.get_atoms())
                    atoms_name = [x.get_name() for x in cycle]

                    try:
                        cycle = cycle[atoms_name.index("CG"):]
                        if residue.get_resname() == "PRO":
                            cycle = [elem for elem in cycle if (elem.get_name() not in ["C", "O"])]
                    except ValueError:
                        pass

                    del atoms_name
                    cihbsList.append(cycle)

                try:
                    self.setNeighbourSearch(residue)
                    # проверяем на наличие простых CN, CS, CO и т.д....
                    if residue.get_resname() in BaseEnums.Groups.singleResidues.value.keys():
                        try:
                            atom = BaseEnums.Groups.singleResidues.value.get(residue.get_resname())
                            atom_obj = list(residue.get_atoms())[
                                [element.get_name() for element in residue.get_atoms()].index(atom)]
                            neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)

                            cihbsList.append(neighbour)
                        except ValueError:
                            pass

                    # проверяем на наличие триплетов OCO, NCN и т.д....
                    if residue.get_resname() in BaseEnums.Groups.complexResidues.value.keys():
                        try:
                            atom = BaseEnums.Groups.complexResidues.value.get(residue.get_resname())
                            atom_obj = list(residue.get_atoms())[
                                [element.get_name() for element in residue.get_atoms()].index(atom)]
                            neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)
                            neighbour = neighbour[[element.get_name() for element in neighbour].index(atom):]

                            cihbsList.append(neighbour)
                        except ValueError:
                            pass
                except IndexError:
                    pass

        # добавляем триплеты NCO из основной цепи
        tmp = [i for i in self.mainChain if (i.get_name() != "CA")]
        cihbsList.append(tmp)

        # созраняем в словарь для удобной и быстрой работы алгоритма
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

    # получаем номер аминокислоты
    def getResseq(self, residue):
        resseq = [st.split('=') for st in str(residue).split(' ')]
        return resseq[4][1]

    # посик ближайших атомов в ССИВС
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

    # соединение атомов основной цепи при помощи физических операторов
    def connectPhysicalOperators(self):
        self.setNeighbourSearch(self.innerCIHBS)

        # получаем пентофрагмент
        for i in range(4, len(self.alphaCarbonChain)):
            group_0 = self.alphaCarbonChain[i].get_parent()
            group_4 = self.alphaCarbonChain[i - 4].get_parent()

            # получаем внутренние ССИВС для i-го элемента

            groupCIHBS = self.innerCIHBSDict.get(group_0)

            def getTargetElementInResidue(group, target):
                return list(group.get_atoms())[list(i.get_name() for i in group.get_atoms()).index(target)]

            i_0_element = getTargetElementInResidue(group_0, "N")
            i_4_element = getTargetElementInResidue(group_4, "O")

            groupCIHBS = groupCIHBS[:groupCIHBS.index(i_0_element)]

            # иначе соединяем 0 и 4
            if math.dist(i_4_element.get_coord(), i_0_element.get_coord()) <= 4:
                self.connectAtoms(i_0_element, i_4_element, BaseEnums.CHIBSBond.physicalOperator)

            # соединяем с остатком
            closest = self.findClosest(groupCIHBS, i_4_element)

            if groupCIHBS and closest and math.dist(closest.get_vector(), i_4_element.get_vector()) <= 4 \
                    and group_4.get_resname() not in BaseEnums.Groups.antiConnectors.value:
                self.connectAtoms(closest, i_4_element, BaseEnums.CHIBSBond.acceptorAcceptor
                if (i_4_element.element == closest.element)
                else BaseEnums.CHIBSBond.acceptorAcceptor)

    # удаляем из группы атомы-ионы
    def checkForIons(self, group, Ion):
        if (Ion == BaseEnums.Atoms.donor):
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.plusIONResidue.value.keys()
                     and i.get_parent() not in BaseEnums.Groups.plusIONResidue.value.values()]
        else:
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.minusIONResidue.value.keys()
                     and i.get_parent() not in BaseEnums.Groups.minusIONResidue.value.values()]
        return group

    # проверяем связь на условие ретроспективности
    def checkRetrospective(self, residue_resseq, neighbours):
        return [i for i in neighbours if int(self.getResseq(i.get_parent())) < int(residue_resseq)]

    # соединение ССИВС между аминокислотами
    def connectResidueCIHBS(self):

        cleanedCIHBS = [i for i in self.innerCIHBS
                        if i.element != "C"
                        and i.get_name() != "N"]

        self.setNeighbourSearch(cleanedCIHBS)
        for elem in cleanedCIHBS:

            neighbours = self.findTargetNeigbourExept(elem, 3.7, ["N", "O", "S"])

            # проверяем на ионы в группе
            neighbours = self.checkForIons(neighbours, BaseEnums.Atoms.donor
            if (elem.element == BaseEnums.Atoms.donor.value) else BaseEnums.Atoms.acceptor)

            # проверяем на ретроспективность ССИВС
            neighbours = self.checkRetrospective(self.getResseq(elem.get_parent()), neighbours)

            if (neighbours != []):
                for neighbour in neighbours:
                    self.connectAtoms(neighbour, elem,
                                      BaseEnums.CHIBSBond.acceptorAcceptor
                                      if (elem.element == neighbour.element and elem.element == "O")
                                      else BaseEnums.CHIBSBond.acceptorAcceptor)
