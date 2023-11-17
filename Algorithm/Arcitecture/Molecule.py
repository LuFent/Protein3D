class Molecule:
    def __init__(self, atoms, identifier):
        self._atoms = atoms
        self._identifier = identifier

    def get_atom_coords(self):
        return [atom.Get_Coords() for atom in self._atoms]

    def get_atoms(self):
        return self._atoms

    def get_length(self):
        return len(self._atoms)

    def add_atom(self, atom):
        self._atoms.append(atom)

    def set_atom(self, atom_arr):
        self._atoms = atom_arr
