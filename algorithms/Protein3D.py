from .Arcitecture import CIHBS, CIHBS_selector
from .CustomPDBParser import PDBParser
from Bio.PDB import Selection, PDBIO, Select
import pprint

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"


class StructureVisualisation:
    carbon_bond_type_number = 4
    PDBParser = PDBParser(QUIET=True)
    io = PDBIO()

    def __init__(self, structure_id, structure_file):

        self.structure = self.PDBParser.get_structure(structure_id, structure_file)

        self.cihbs_obj = CIHBS.BaseCIHBS()
        self.inner_cihbs_obj = CIHBS.InnerCIHBS()
        self.outer_cihbs_obj = CIHBS.OuterCIHBS()
        self.physical_cihbs_obj = CIHBS.PhysicalOperators()
        self.cihbs_selector_obj = CIHBS_selector.SelectorCIHBS()

        self.cihbs_obj.setParameters(structure=self.structure)

        self.atoms_ids = [atom.get_serial_number() for atom in
                          Selection.unfold_entities(self.structure, 'A')]

        self.extra_bonds = []
        self.hide_atoms = []
        self.show_atoms = []

    @classmethod
    def get_structure(cls, structure_id, structure_file):
        return cls.PDBParser.get_structure(structure_id, structure_file)

    @classmethod
    def clean_file(cls, structure, file):
        cls.io.set_structure(structure)
        cls.io.save(file, select=NotDisordered())

    def flush_mask(self):
        self.hide_atoms = []
        self.extra_bonds = []
        self.show_atoms = []

    def get_mask(self):
        return {
            "hideAtoms": self.hide_atoms,
            "showAtoms": self.show_atoms,
            "extraBonds": self.extra_bonds,
        }

    def get_alpha_carbon_chain(self):
        self.flush_mask()

        show_ids = []
        extra_bonds = []

        previous_atom = None

        for model in self.structure:
            # Iterate through each chain in the model
            for chain in model:
                # Iterate through each residue in the chain
                for residue in chain:

                    # Iterate through each atom in the residue
                    for atom in residue:

                        if atom.get_name() == "CA":
                            show_ids.append(atom.get_serial_number())
                            if previous_atom is None:
                                previous_atom = atom.get_serial_number()
                            else:
                                extra_bonds.append(
                                    (previous_atom, atom.get_serial_number(), self.carbon_bond_type_number))
                                previous_atom = atom.get_serial_number()
                previous_atom = None

        self.show_atoms = show_ids
        self.hide_atoms = None
        self.extra_bonds = {"alpha_carbon_skeleton": {
            "alpha_carbon_skeleton": extra_bonds
        }
        }

    def getCIHBS(self):
        self.flush_mask()

        self.inner_cihbs_obj.checkInnerGroups(self.structure, self.cihbs_obj.getMainChain())

        self.physical_cihbs_obj.setParameters(args=(self.inner_cihbs_obj.getInnerCIHBSDict(),
                                                    self.inner_cihbs_obj.getInnerCIHBS(),
                                                    self.cihbs_obj.getCarbonChain()))
        self.physical_cihbs_obj.connectPhysicalOperators()
        new_cihbs = self.physical_cihbs_obj.getNewCIHBS()

        self.outer_cihbs_obj.setParameters(self.inner_cihbs_obj.getInnerCIHBS())
        self.outer_cihbs_obj.connectResidueCIHBS()
        new_cihbs.extend(self.outer_cihbs_obj.getNewCIHBS())

        self.cihbs_selector_obj.setParameters(new_cihbs)
        self.extra_bonds = self.cihbs_selector_obj.get_selector_map()
        self.hide_atoms = None
        self.show_atoms = self.inner_cihbs_obj.getInnerCIHBS_serial_number()

    def test_alg(self):
        """self.flush_mask()
        atom_mask = self.cihbsObj.getPropertiesColoredAtoms(self.structure)"""
        pass

    def execute_algorithm(self, alg):
        if alg == "alpha_carbon_skeleton":
            self.get_alpha_carbon_chain()
            return self.get_mask()

        elif alg == "CIHBS":
            self.getCIHBS()
            return self.get_mask()

        elif alg == "test_alg":
            self.test_alg()
            return self.get_mask()

        else:
            return None
