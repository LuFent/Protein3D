
from .Arcitecture import CIHBS
from .CustomPDBParser import PDBParser
from Bio import PDB as pdb
from Bio.PDB import Selection

from itertools import chain, count


class StructureVisualisation:
    cihbsObj = CIHBS.CIHBS()
    PDBParser = PDBParser(QUIET=True)

    def __init__(self, structure_id, structure_file):

        structure = self.PDBParser.get_structure(structure_id, structure_file)
        self.structure = structure
        self.atoms_ids = [atom.get_serial_number() for atom in
               Selection.unfold_entities(self.structure, 'A')]

        self.extra_bonds = []
        self.hide_atoms = []

    def flush_mask(self):
        self.hide_atoms = []
        self.extra_bonds = []

    def get_mask(self):
        return {
            "hide_atoms": self.hide_atoms,
            "extra_bonds": self.extra_bonds
        }

    def get_alpha_carbon_chain(self):
        self.flush_mask()

        hide_ids = [atom.get_serial_number() for atom in
               Selection.unfold_entities(self.structure, 'A') if atom.get_name() != "CA"]

        self.hide_atoms = hide_ids

    def getCHIBS(self):
        self.flush_mask()
        innerCIHBS = self.cihbsObj.checkInnerGroups(self.structure)
        CIHBS_ids = set(chain.from_iterable(innerCIHBS))
        # self.cihbsObj.checkMultiGroups()
        self.hide_atoms = [atom for atom in self.atoms_ids if atom not in CIHBS_ids]