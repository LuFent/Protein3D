
class Molecule():
    def __init__(self, atoms, identifier):
        self._atoms = atoms
        self._identifier = identifier

    def Get_Atom_Coords(self):
        return [atom.Get_Coords() for atom in self._atoms]

    def Get_Atoms(self):
        return self._atoms

    def Get_Length(self):
        return len(self._atoms)

    def Add_Atom(self, atom):
        self._atoms.append(atom)

    def Set_Atom(self, atom_arr):
        self._atoms = atom_arr