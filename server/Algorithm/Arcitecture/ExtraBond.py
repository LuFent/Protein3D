from .BaseEnums import CHIBSBond

from pprint import pprint


class ExtraBonds:
    def __init__(self, atom_1, atom_2, bond_type, bond_energy):
        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.bond_type = bond_type
        self.bond_energy = bond_energy

    def getAtom1(self):
        return self.atom_1

    def getAtom2(self):
        return self.atom_2

    def getBondType(self):
        return self.bond_type

    def getBondEnergy(self):
        return self.bond_energy

    def getAll(self):
        return [self.atom_1.get_serial_number(), self.atom_2.get_serial_number(), self.bond_type, self.bond_energy]

    def getAll_no_Bondtype(self):
        return [self.atom_1.get_serial_number(), self.atom_2.get_serial_number(), self.bond_energy]


class ExtraBondStorage:
    def __init__(self):
        self.bonds = []

    def getBonds(self):
        return self.bonds

    def addBond(self, bond):
        self.bonds.append(bond)

    def extendBond(self, bonds):
        self.bonds.extend(bonds)

    def mapping_NH_O(self, distance, bond_type):
        return [bond.getAll_no_Bondtype() for bond in self.bonds if
                abs(int(bond.getAtom1().get_parent().get_id()[1]) - int(
                    bond.getAtom2().get_parent().get_id()[1])) == distance
                and bond.getAtom1().get_name() == 'N'
                and bond.getAtom2().get_name() == 'O'
                and bond.getBondType() == bond_type]

    def mapping_R_O(self, distance, bond_type):
        return [bond.getAll_no_Bondtype() for bond in self.bonds if
                abs(int(bond.getAtom1().get_parent().get_id()[1])) - int(
                    bond.getAtom2().get_parent().get_id()[1]) == distance
                and len(bond.getAtom1().get_name()) >= 2
                and bond.getAtom2().get_name() == 'O'
                and bond.getBondType() == bond_type]

    def mapping_R_R(self, distance, bond_type):
        return [bond.getAll_no_Bondtype() for bond in self.bonds if
                abs(int(bond.getAtom1().get_parent().get_id()[1])) - int(
                    bond.getAtom2().get_parent().get_id()[1]) == distance
                and len(bond.getAtom1().get_name()) >= 2
                and len(bond.getAtom2().get_name()) >= 2
                and bond.getBondType() == bond_type]

    # получение мостов между аминокислотами
    def getBridges(self):
        return {'Polypeptide': [bond.getAll_no_Bondtype() for bond in self.bonds
                                if bond.getBondType() == CHIBSBond.residueBridge.value.get("type")]}

    # получение остальных ССИВС
    def getOthers(self, bond_type):
        return [bond.getAll_no_Bondtype() for bond in self.bonds
                if bond.getBondType() != CHIBSBond.physicalOperator.value.get("type")
                and bond.getBondType() == bond_type
                and ((len(bond.getAtom1().get_name()) >= 2 and bond.getAtom2().get_name() != 'O')
                     or (len(bond.getAtom1().get_name()) == 1 and len(bond.getAtom2().get_name()) >= 2))]

    def serialize(self):
        if self.bonds:
            bonds = {}

            bonds['HydrogenAA'] = {
                'Others': [i for i in self.getOthers(CHIBSBond.acceptorAcceptor.value.get("type"))],

                # получение ССИВС между iм остатком и i-n кислородом пентофрагментов
                'R_Oi_{}'.format(3): [i for i in self.mapping_R_O(3, CHIBSBond.acceptorAcceptor.value.get("type"))],

                # получение ССИВС между iм остатком и i-n остатком
                'R_Ri_{}'.format(3): [i for i in self.mapping_R_R(3, CHIBSBond.acceptorAcceptor.value.get("type"))]}
            bonds['HydrogenDA'] = {
                'Others': [i for i in self.getOthers(CHIBSBond.donorAcceptor.value.get("type"))],
                'R_Oi_{}'.format(3): [i for i in self.mapping_R_O(3, CHIBSBond.donorAcceptor.value.get("type"))],
                'R_Ri_{}'.format(3): [i for i in self.mapping_R_R(3, CHIBSBond.donorAcceptor.value.get("type"))]}

            # получение физических операторов между i и i-n элементами пентофрагментов
            bonds['Physics'] = {'NH_Oi_{}'.format(3): [i for i in self.mapping_NH_O(3, CHIBSBond.physicalOperator.value.get("type"))]}

            for cycle in range(4, 7):
                bonds['HydrogenAA'].update({
                    'R_Oi_{}'.format(cycle): [i for i in self.mapping_R_O(cycle, CHIBSBond.acceptorAcceptor.value.get("type"))],
                    'R_Ri_{}'.format(cycle): [i for i in self.mapping_R_R(cycle, CHIBSBond.acceptorAcceptor.value.get("type"))]})
                bonds['HydrogenDA'].update({
                    'R_Oi_{}'.format(cycle): [i for i in self.mapping_R_O(cycle, CHIBSBond.donorAcceptor.value.get("type"))],
                    'R_Ri_{}'.format(cycle): [i for i in self.mapping_R_R(cycle, CHIBSBond.donorAcceptor.value.get("type"))]})
                bonds['Physics'].update({'NH_Oi_{}'.format(cycle): [i for i in self.mapping_NH_O(cycle, CHIBSBond.physicalOperator.value.get("type"))]})

            bonds['Bridges'] = self.getBridges()
            # pprint(bonds)
            return bonds
        raise ValueError("Возвращается пустой объект.")
