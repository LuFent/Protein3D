import math

from Bio.PDB import Selection, Atom, NeighborSearch

from Arcitecture import BaseEnums

class CIHBS:
    def __init__(self):
        self.newCIHBS = []

    def setParameters(self, structure):
        self.alphaCarbonChain = [atom for atom in Selection.unfold_entities(structure, 'A') if atom.get_name() == "CA"]
        self.mainChain = [atom for atom in Selection.unfold_entities(structure, 'A') if
                          ((atom.get_name() == "CA" or len(atom.get_name()) == 1)
                           and atom.get_parent().get_resname() != "HOH")]
        self.ns = NeighborSearch(Selection.unfold_entities(structure, 'A'))

        self.nsStructure = NeighborSearch(Selection.unfold_entities(structure, 'A'))

    def getMainChain(self):
        return self.mainChain

    def getCarbonChain(self):
        return self.alphaCarbonChain

    def getInnerCIHBS(self):
        return sum(self.innerCIHBS, [])

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
            if residue.get_resname() != "HOH":
                #checking for benzene rings and adding them
                if residue.get_resname() in BaseEnums.Groups.sycleResidues.value:
                    cycle = list(residue.get_atoms())
                    atomsName = [x.get_name() for x in cycle]

                    try:
                        cycle = cycle[atomsName.index("CG"):]
                        if residue.get_resname() == "PRO":
                            cycle = [elem for elem in cycle if (elem.get_name() not in ["C", "O"])]
                    except ValueError:
                        #print("Не найден опорный атом")
                        pass

                    #print("atoms cycle", cycle)
                    del atomsName
                    cihbsList.append(cycle)

                self.setNeighbourSearch(residue)
                #simple residues CN, CS, CO and etc...
                if residue.get_resname() in BaseEnums.Groups.singleResidues.value.keys():
                    try:
                        atom = BaseEnums.Groups.singleResidues.value.get(residue.get_resname())
                        atom_obj = list(residue.get_atoms())[ [element.get_name() for element in residue.get_atoms()].index(atom) ]
                        neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)


                        #print("Neighbours simple group = ", len(neighbour), neighbour)
                        cihbsList.append(neighbour)
                    except ValueError:
                        #print("Не найден опорный атом", BaseEnums.Groups.singleResidues.value.get(residue.get_resname()))
                        pass

                #all troplets OCO, NCN and etc...
                if residue.get_resname() in BaseEnums.Groups.complexResidues.value.keys():
                    try:
                        atom = BaseEnums.Groups.complexResidues.value.get(residue.get_resname())
                        atom_obj = list(residue.get_atoms())[ [element.get_name() for element in residue.get_atoms()].index(atom) ]
                        neighbour = self.findTargetNeigbour(atom_obj, 2.0, BaseEnums.Atoms.allAtomsElements.value)
                        neighbour = neighbour[ [element.get_name() for element in neighbour].index(atom): ]

                        #print("Neighbours complex group = ", len(neighbour), neighbour)
                        cihbsList.append(neighbour)
                    except ValueError:
                        #print("Не найден опорный атом", BaseEnums.Groups.complexResidues.value.get(residue.get_resname()))
                        pass

        #NCO
        tmp = [i for i in self.mainChain if(i.get_name() != "CA")]
        cihbsList.append(tmp)

        self.innerCIHBS = cihbsList



    def removeCarbon(self, group):
        tmp = []
        for elem in group:
            #print("elem.element in removeCarbon", elem.element)
            if elem.get_name()[0] != "C":
                tmp.append(elem)
        return tmp

    def getResidueAtoms(self, atom):
        #print("parent = ", atom.get_parent())
        return list(atom.get_parent().get_atoms())

    def getResseq(self, residue):
        resseq = [st.split('=') for st in str(residue).split(' ')]
        return resseq[4][1]

    def findInnerCIHBS(self, residue):
        atoms = [i for i in self.getInnerCIHBS() if i.get_parent() == residue]
        return atoms

    def findClosest(self, atoms, target):
        atoms = [i for i in atoms if(i.element != "C"
                                     and i.get_parent() not in BaseEnums.Groups.doubleBondAcid.value.keys()
                                     and i.get_name() not in BaseEnums.Groups.doubleBondAcid.value.values()
                                     and  not i.is_disordered() )]
        try:
            distances = [math.dist(i.get_coord(), target.get_coord()) for i in atoms]
            return atoms[distances.index(min(distances))]
        except ValueError:
            return None

    def findAlternativePhysicalOperator(self, target_element, find_element):
        self.setNeighbourSearch(sum(self.innerCIHBS, []))
        neighbours = [i for i in self.findTargetNeigbour(target_element, 4.0, find_element)
                      if len(i.get_name()) > 1]
        print("chain connection failed", str(target_element), target_element.get_parent())
        if neighbours != []:
            print("no N found, found alternative variant",
                  [i.get_parent() for i in neighbours], neighbours)
            print("closest is", self.findClosest(neighbours, target_element).get_parent(),
                  self.findClosest(neighbours, target_element))
            self.connectAtoms(self.findClosest(neighbours, target_element), target_element,
                              BaseEnums.CHIBSBond.physicalOperator)

    def connectPhysicalOperators(self):

        # получаем пентофрагмент
        for i in range(5, len(self.alphaCarbonChain)):
            group = self.alphaCarbonChain[i-5:i]

            if group[4].get_parent().get_resname() not in BaseEnums.Groups.antiConnectors.value:
                # получаем внутренние ССИВС для i-го элемента
                #print("for alpha group", [j.get_parent() for j in group])
                groupCIHBS = self.findInnerCIHBS(group[4].get_parent())
                #print("groupCIHBS", [i.get_parent() for i in groupCIHBS], groupCIHBS)

                i_0_element = groupCIHBS[ [i.get_name() for i in groupCIHBS].index("N") ]
                i_4_element = self.findInnerCIHBS(group[0].get_parent())
                i_4_element = i_4_element[ [i.get_name() for i in i_4_element].index("O") ]

                groupCIHBS = groupCIHBS[:groupCIHBS.index(i_0_element)]

                #если SER, THR перекрывают, то соединяем остаток и 3
                #if(groupCIHBS and len(groupCIHBS) == 1 and math.dist(groupCIHBS[0].get_coord(), i_0_element.get_coord()) <= 4):
                #    print("connecting in res SER THR = ", math.dist(groupCIHBS[0].get_coord(), i_0_element.get_vector()), groupCIHBS[0].get_parent() )
                #    self.connectAtoms(i_3_element, groupCIHBS[0], BaseEnums.CHIBSBond.outerCIHBS)

                #иначе соединяем 0 и 4
                if math.dist(i_4_element.get_coord(), i_0_element.get_coord()) <= 4:
                    #print("connecting in chain = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))
                    self.connectAtoms(i_0_element, i_4_element, BaseEnums.CHIBSBond.physicalOperator)
                #else:
                    #self.findAlternativePhysicalOperator(i_4_element, 'N')
                    #self.findAlternativePhysicalOperator(i_0_element, 'O')

                # соединяем с остатком
                # сделать проверку на двойную связь

                if groupCIHBS and self.findClosest(groupCIHBS, i_4_element) and math.dist(self.findClosest(groupCIHBS, i_4_element).get_vector(), i_4_element.get_vector()) <= 4:
                    #print("connecting residual_4 to ", [i.get_parent() for i in groupCIHBS], "with DIST = ", math.dist(i_0_element.get_vector(), i_4_element.get_vector()))
                    self.connectAtoms(self.findClosest(groupCIHBS, i_4_element), i_4_element, BaseEnums.CHIBSBond.physicalOperator)

                #print("\n")
        #print(self.newCIHBS)
        return self.newCIHBS



    def removeNitrogen(self, group):
        tmp = []
        for elem in group:
            #print("elem.element in removeCarbon", elem.element)
            if( (elem.get_name()[0] == "N" and len(elem.get_name()) != 1) or elem.get_name()[0] != "N"):
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
        return [i for i in group if( i not in self.mainChain)]

    def sordByDistance(self, target, group):
        pass

    def connectResidueCIHBS(self):

        #чистим
        cleanedCIHBS = [i for i in self.getInnerCIHBS()
                         if i.element != "C"
                         and i.get_name() != "N"]
        #print(cleanedCIHBS)

        gotInitiatorAcid = False
        for elem in cleanedCIHBS:
            #print("for res", [i.get_parent() for i in group], group)

            #if(elem.get_parent().get_resname() in BaseEnums.Groups.CIHBSInitialisators.value):
            #    gotInitiatorAcid = True

            self.setNeighbourSearch( cleanedCIHBS ) #sum(map, [])
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
