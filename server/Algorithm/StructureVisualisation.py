from .Arcitecture import CIHBS, CIHBS_selector
from .CustomPDBParser import PDBParser
from Bio.PDB import Selection, PDBIO, Select
from .Algorithm import Algorithm, AlgorithmsStorage
from .Mask import Mask

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"


class StructureVisualisation:
    carbon_bond_type_number = 4
    PDBParser = PDBParser(QUIET=True)
    io = PDBIO()


    def __init__(self, structure_id, structure_file, algorithms_storage, filename):
        self.structure = self.PDBParser.get_structure(structure_id, structure_file)
        self.mask = Mask()
        self.algorithms_storage = algorithms_storage
        self.filename = filename
        self.atoms_ids = [atom.get_serial_number() for atom in
                          Selection.unfold_entities(self.structure, 'A')]


    @classmethod
    def get_structure(cls, structure_id, structure_file):
        return cls.PDBParser.get_structure(structure_id, structure_file)

    @classmethod
    def clean_file(cls, structure, file):
        cls.io.set_structure(structure)
        cls.io.save(file, select=NotDisordered())


    def flush_mask(self):
        self.mask.flush()

    def get_mask(self):
        return self.mask.serialize()


    def execute_algorithm(self, alg_code_name):
        alg = self.algorithms_storage[alg_code_name]
        self.flush_mask()
        return alg.execute(self.structure)


    def __str__(self):
        return f"Structure -- {self.filename}"