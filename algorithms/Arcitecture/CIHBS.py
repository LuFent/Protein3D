from Bio.PDB import Selection, Atom, NeighborSearch

from Arcitecture import BaseEnums


class CIHBS:
    def __init__(self):
        self._atomO = Atom.Atom("O", None, None, None, '', "O", None, "O")
        self._atomCG = Atom.Atom("CG", None, None, None, '', "CG", None, "C")

    def findNeighbours(self, residue, atom, condotion):
        # поиск соседей с опорным атомом включительно
        ns = NeighborSearch(Selection.unfold_entities(residue, 'A'))

        neighbours = [elem for elem in ns.search(atom.get_coord(), 2.0, 'A')]
        if condotion:
            neighbours = [elem for elem in ns.search(atom.get_coord(), 2.0, 'A') if (elem.get_name()[0] != condotion or elem == atom)]
        return neighbours

    def checkInnerGroups(self, structure):

        cihbsList = []
        for residue in structure.get_residues():
            #print("\nfor residue = ", residue.get_resname())

            # проверка на бензольные кольца и их добавление
            if (residue.get_resname() in BaseEnums.Groups.sycleResiduces.value):
                #print("проверка на бензольные кольца и их добавление")
                cycle = list(residue.get_atoms())
                atomsName = [x.get_name() for x in cycle]

                try:
                    cycle = cycle[atomsName.index(self._atomCG.get_name()):]
                    cycle = [atom.get_serial_number() for atom in cycle] #!!!!!!!

                except ValueError:
                    print("Не найден опорный атом CG")
                    pass
                #print("atoms cycle", cycle)

                del atomsName
                cihbsList.append(tuple(cycle))

            for atom in residue:
                #print("for atom ", atom.get_name(), type(atom.get_name()))

                # проверка простых групп
                if (residue.get_resname() in BaseEnums.Groups.singleResidues.value):
                    #print("проверка простых групп")
                    if (len(atom.get_name()) > 1 and atom.get_name()[0] in BaseEnums.Atoms.singleTargetAtoms.value):
                        neighbour = self.findNeighbours(residue, atom, None)
                        #print("Neighbours = ", len(neighbour), neighbour)

                        if (len(neighbour) > 2):
                            neighbour.pop()
                            #print("Превышение количества соседей для ССИВС в простой группе")

                        # добавляем
                        neighbour = [atom.get_serial_number() for atom in neighbour] #!!!!!!!!!
                        #print("Neighbours after cut = ", len(neighbour), neighbour)
                        cihbsList.append(tuple(neighbour))

                # проверка резонансных групп, проверка связей N-C-O, O-C-O, N-C-N и тупиковые C-O
                if (atom.element == "C"):
                    #print("проверка резонансных групп")
                    neighbour = self.findNeighbours(Selection.unfold_entities(structure, 'A'), atom, "C")

                    if (len(neighbour) > 3 or self._atomO.element in [elem.element for elem in neighbour]):
                        #print("Neighbours resonance = ", len(neighbour), neighbour)

                        # добавляем
                        neighbour = [atom.get_serial_number() for atom in neighbour] #!!!!!!!!!
                        cihbsList.append(tuple(neighbour))

        return cihbsList

    def checkMultiGroups(self):
        """
            ALA, VAL, LEU, ILE, PHE не образуют

            В зависимости от возможной молекулярной функции в составе ССИВС активные
            элементы были разделены на ряд групп [1–4] (табл. 6.3).
                1. Инициаторы ССИВС: пролин (Pro), метионин (Met).
                2. Молекулярные клапаны: (Ser), (Thr), (Asp), (Glu), (Asn), (Glu), (Lys), (Arg).
                3. Элементы задержки сигнала: (His), (Tyr), (Trp).
                4. Элементы инверсии сигнала: 3-метилгистидин (3-me-His), (Tyr), (Trp).
        """
        pass