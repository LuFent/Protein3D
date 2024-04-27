from server import run_server
from algorithms.Algorithm import Algorithm, AlgorithmsStorage
from algorithms.Mask import Mask
from algorithms.Arcitecture import CIHBS, CIHBS_selector


class AlphaCarbonChainAlgorithm(Algorithm):
    code_name = "alpha_carbon_skeleton"
    label = "Alpha-Carbon Skeleton"


    def _execute(self, structure) -> Mask:
        carbon_bond_type_number = 4
        mask = Mask()
        show_ids = []
        extra_bonds = []

        previous_atom = None

        for model in structure:
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
                                    (previous_atom, atom.get_serial_number(), carbon_bond_type_number))
                                previous_atom = atom.get_serial_number()
                previous_atom = None

        mask.show_atoms = show_ids
        mask.hide_atoms = None
        mask.extra_bonds = {"alpha_carbon_skeleton": {
            "alpha_carbon_skeleton": extra_bonds
            }
        }
        return mask


class CIHBSAlgorithm(Algorithm):
    code_name = "CIHBS"
    label = "Ion-Hydrogen Bonds System"

    def __init__(self):
        self.cihbs_obj = CIHBS.BaseCIHBS()
        self.inner_cihbs_obj = CIHBS.InnerCIHBS()
        self.outer_cihbs_obj = CIHBS.OuterCIHBS()
        self.physical_cihbs_obj = CIHBS.PhysicalOperators()
        self.cihbs_selector_obj = CIHBS_selector.SelectorCIHBS()


    def _execute(self, structure) -> Mask:
        mask = Mask()
        self.cihbs_obj.setParameters(structure=structure)
        self.inner_cihbs_obj.checkInnerGroups(structure, self.cihbs_obj.getMainChain())

        self.physical_cihbs_obj.setParameters(args=(self.inner_cihbs_obj.getInnerCIHBSDict(),
                                                    self.inner_cihbs_obj.getInnerCIHBS(),
                                                    self.cihbs_obj.getCarbonChain()))
        self.physical_cihbs_obj.connectPhysicalOperators()
        new_cihbs = self.physical_cihbs_obj.getNewCIHBS()

        self.outer_cihbs_obj.setParameters(self.inner_cihbs_obj.getInnerCIHBS())
        self.outer_cihbs_obj.connectResidueCIHBS()
        new_cihbs.extend(self.outer_cihbs_obj.getNewCIHBS())

        self.cihbs_selector_obj.setParameters(new_cihbs)
        mask.extra_bonds = self.cihbs_selector_obj.get_selector_map()
        mask.hide_atoms = None
        mask.show_atoms = self.inner_cihbs_obj.getInnerCIHBS_serial_number()
        return mask


run_server(AlphaCarbonChainAlgorithm, CIHBSAlgorithm)