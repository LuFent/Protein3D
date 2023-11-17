class Atom():

    def __init__(self, hparams):
        self.data_col = hparams[0]
        self._atom_name = hparams[1]
        self._residue = hparams[2]
        self._identifier = hparams[3]
        self._sequence_number = hparams[4]
        self._x = hparams[5]
        self._y = hparams[6]
        self._z = hparams[7]
        self._occupancy = hparams[8]
        self._b_factor = hparams[9]
        self._elementsymbol = hparams[10]

    def Get_atom_name(self):
        return self._atom_name

    def Get_identifier(self):
        return self._identifier

    def Get_Coords(self):
        return self._x, self._y, self._z

    def Create_Bond(self, atom_connect):
        pass



"""
    def __init__(self, atom_name, residue , identifier, sequence_number,
                 x, y, z, occupancy, b_factor, elementsymbol):

        self._atom_name = atom_name
        self._residue = residue
        self._identifier = identifier
        self._sequence_number = sequence_number
        self._x = x
        self._y = y
        self._z = z
        self._occupancy = occupancy
        self._b_factor = b_factor
        self._elementsymbol = elementsymbol

        self.Connected_atoms = {}
    """
