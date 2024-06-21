import math

from Bio.PDB import Selection, NeighborSearch

from . import BaseEnums


class BaseCIHBS:
    def __init__(self):
        self.new_cihbs = []

    def setParameters(self, structure):
        self.alpha_carbon_chain = [atom for atom in Selection.unfold_entities(structure, 'A')
                                   if atom.get_name() == "CA"]
        """self.main_chain_dict = {atom.get_parent(): atom.get_name() for atom in Selection.unfold_entities(structure, 'A')
                              if ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                                  and atom.get_parent().get_resname() != "HOH")}"""

        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))
        self.setMainChain(structure)

    def getNewCIHBS(self):
        return self.new_cihbs

    def setMainChain(self, structure):
        self.main_chain = [atom for atom in Selection.unfold_entities(structure, 'A')
                           if ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                               and atom.get_parent().get_resname() != "HOH")]

    # получить номера атомов основной цепи
    def getMainChain_serial_number(self):
        return [i.get_serial_number() for i in self.main_chain]

    # получить атомы основной цепи
    def getMainChain(self):
        return self.main_chain

    # получить альфауглерод
    def getCarbonChain(self):
        return self.alpha_carbon_chain


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
        self.new_cihbs.append([atom1, atom2, bond])


class InnerCIHBS(BaseCIHBS):
    def __init__(self):
        super().__init__()
        self.inner_cihbs_dict = {}

    # получить номера атомов простых ССИВС
    def getInnerCIHBS_serial_number(self):
        return [atom.get_serial_number() for atom in self.inner_cihbs]

    # получить номера атомов простых ССИВС
    def getInnerCIHBS(self):
        return self.inner_cihbs

    # получить словарь простых ССИВС
    def getInnerCIHBSDict(self):
        return self.inner_cihbs_dict

    def checkInnerGroups(self, structure, main_chain):

        cihbs_list = []
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
                    cihbs_list.append(cycle)

                try:
                    self.setNeighbourSearch(residue)
                    # проверяем на наличие простых CN, CS, CO и т.д....
                    if residue.get_resname() in BaseEnums.Groups.singleResidues.value.keys():
                        try:
                            atom = BaseEnums.Groups.singleResidues.value.get(residue.get_resname())
                            atom_obj = list(residue.get_atoms())[
                                [element.get_name() for element in residue.get_atoms()].index(atom)]
                            neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)

                            cihbs_list.append(neighbour)
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

                            cihbs_list.append(neighbour)
                        except ValueError:
                            pass
                except IndexError:
                    pass

        # добавляем триплеты NCO из основной цепи
        tmp = [i for i in main_chain if (i.get_name() != "CA")]
        cihbs_list.append(tmp)

        # созраняем в словарь для удобной и быстрой работы алгоритма
        self.inner_cihbs = sum(cihbs_list, [])
        for i in self.inner_cihbs:
            if i.get_parent() in self.inner_cihbs_dict.keys():
                tmp_list = self.inner_cihbs_dict.get(i.get_parent())
                tmp_list.append(i)
                self.inner_cihbs_dict.update({i.get_parent(): tmp_list})
                del tmp_list
            else:
                self.inner_cihbs_dict[i.get_parent()] = [i]


class PhysicalOperators(BaseCIHBS):
    def __init__(self):
        super().__init__()

    def setParameters(self, args):
        self.inner_cihbs_dict, self.inner_cihbs, self.alpha_carbon_chain = args

    @staticmethod
    # посик ближайших атомов в ССИВС
    def findClosest(atoms, target):
        atoms = [i for i in atoms if (i.element != "C"
                                      and i.get_parent() not in BaseEnums.Groups.doubleBondAcid.value.keys()
                                      and i.get_name() not in BaseEnums.Groups.doubleBondAcid.value.values()
                                      and not i.is_disordered())]
        try:
            distances = [math.dist(i.get_coord(), target.get_coord()) for i in atoms]
            return atoms[distances.index(min(distances))]
        except ValueError:
            return None

    def checkResidue(self, group_cihbs, i_4_element, group_4_resname):
        closest = self.findClosest(group_cihbs, i_4_element)

        if group_cihbs and closest and math.dist(closest.get_vector(), i_4_element.get_vector()) <= 4 \
                and group_4_resname not in BaseEnums.Groups.antiConnectors.value:
            self.connectAtoms(closest, i_4_element, BaseEnums.CHIBSBond.acceptorAcceptor
            if (i_4_element.element == closest.element)
            else BaseEnums.CHIBSBond.acceptorAcceptor)

    # соединение атомов основной цепи при помощи физических операторов
    def connectPhysicalOperators(self):
        self.setNeighbourSearch(self.inner_cihbs)

        # получаем пентофрагмент
        for i in range(4, len(self.alpha_carbon_chain)):
            group_0 = self.alpha_carbon_chain[i].get_parent()
            group_4 = self.alpha_carbon_chain[i - 4].get_parent()
            group_3 = self.alpha_carbon_chain[i - 3].get_parent()

            #проверка на наличие аминокислот в пентофрагменте
            if group_0.get_resname() not in BaseEnums.Groups.all_amino_acid_list.value or \
                group_4.get_resname() not in BaseEnums.Groups.all_amino_acid_list.value:
                continue

            diff = int(group_0.get_id()[1]) - int(group_4.get_id()[1])
            #print("for i = ", i, diff, group_0, group_4)
            # проверка на совпадения серийных номеров и возвращение к пентофрагменту
            if diff != 4:
                group_4 = self.alpha_carbon_chain[i - 4 + (diff - 4)].get_parent()
                group_3 = self.alpha_carbon_chain[i - 3 + (diff - 3)].get_parent()
                diff = int(group_0.get_id()[1]) - int(group_4.get_id()[1])
               #print("New diff", diff)
            if diff == 4:
                pass
            else:
                continue

            #print("Checking resseq = ", group_0, group_4)

            # получаем внутренние ССИВС для i-го элемента

            group_cihbs = self.inner_cihbs_dict.get(group_0)

            def getTargetElementInResidue(group, target):
                return list(group.get_atoms())[list(i.get_name() for i in group.get_atoms()).index(target)]

            i_0_element = getTargetElementInResidue(group_0, "N")
            i_4_element = getTargetElementInResidue(group_4, "O")
            i_3_element = getTargetElementInResidue(group_3, "O")

            group_cihbs = group_cihbs[:group_cihbs.index(i_0_element)]

            # соединяем 0 и 4
            if math.dist(i_4_element.get_coord(), i_0_element.get_coord()) <= 4:
                #print(i_0_element.get_parent(), i_4_element.get_parent(), "\n")
                self.connectAtoms(i_0_element, i_4_element, BaseEnums.CHIBSBond.physicalOperator)

            # соединяем с остатком
            self.checkResidue(group_cihbs, i_4_element, group_4.get_resname())

            # оцениваем четырехзвенный цикл для SER и THR
            if group_0.get_resname() in list(BaseEnums.AcidGroups.weaklyPolar.value.values())[0]:
                # грубо заменяем переменную
                i_4_element = i_3_element

                # соединяем с остатком i-3
                self.checkResidue(group_cihbs, i_4_element, group_4.get_resname())


class OuterCIHBS(BaseCIHBS):
    def __init__(self):
        super().__init__()

    def setParameters(self, inner_cihbs):
        self.inner_cihbs = inner_cihbs

    @staticmethod
    # удаляем из группы атомы-ионы
    def checkForIons(group, ion):
        if ion == BaseEnums.Atoms.donor:
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.plusIONResidue.value.keys()
                     and i.get_parent() not in BaseEnums.Groups.plusIONResidue.value.values()]
        else:
            group = [i for i in group if i.get_name() not in BaseEnums.Groups.minusIONResidue.value.keys()
                     and i.get_parent() not in BaseEnums.Groups.minusIONResidue.value.values()]
        return group

    # проверяем связь на условие ретроспективности
    def checkRetrospective(self, residue_resseq, neighbours):
        return [i for i in neighbours if int(i.get_parent().get_id()[1]) < int(residue_resseq)]

    # соединение ССИВС между аминокислотами
    def connectResidueCIHBS(self):

        cleaned_cihbs = [i for i in self.inner_cihbs
                         if i.element != "C"
                         and i.get_name() != "N"]

        self.setNeighbourSearch(cleaned_cihbs)
        for elem in cleaned_cihbs:

            neighbours = self.findTargetNeigbourExept(elem, 3.7, ["N", "O", "S"])

            # проверяем на ионы в группе
            neighbours = self.checkForIons(neighbours, BaseEnums.Atoms.donor
            if (elem.element == BaseEnums.Atoms.donor.value) else BaseEnums.Atoms.acceptor)

            # проверяем на ретроспективность ССИВС
            neighbours = self.checkRetrospective(elem.get_parent().get_id()[1], neighbours)

            if neighbours != []:
                for neighbour in neighbours:
                    self.connectAtoms(neighbour, elem,
                                      BaseEnums.CHIBSBond.acceptorAcceptor
                                      if (elem.element == neighbour.element and elem.element == "O")
                                      else BaseEnums.CHIBSBond.acceptorAcceptor)