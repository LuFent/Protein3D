from server.Algorithm.Arcitecture import CIHBS
from server.Algorithm.Algorithm import Algorithm, Mask


class CIHBSAlgorithm(Algorithm):
    code_name = "CIHBS"
    label = "Ion-Hydrogen Bonds System"
    icon = "CIHBSAlgorithmLogo.svg"

    def __init__(self):
        pass

    def set_cihbs(self):
        self.cihbs_obj = CIHBS.BaseCIHBS()
        self.inner_cihbs_obj = CIHBS.InnerCIHBS()
        self.outer_cihbs_obj = CIHBS.OuterCIHBS()
        self.physical_cihbs_obj = CIHBS.PhysicalOperators()

    def _execute(self, structure) -> Mask:
        self.set_cihbs()
        mask = Mask()
        self.cihbs_obj.setParameters(structure=structure)
        new_cihbs = self.cihbs_obj.getNewCIHBS()

        self.inner_cihbs_obj.checkInnerGroups(structure, self.cihbs_obj.getMainChain())

        self.physical_cihbs_obj.setParameters(args=(self.inner_cihbs_obj.getInnerCIHBSDict(),
                                                    self.inner_cihbs_obj.getInnerCIHBS(),
                                                    self.cihbs_obj.getCarbonChain()))
        self.physical_cihbs_obj.connectPhysicalOperators()
        new_cihbs.extendBond(self.physical_cihbs_obj.getNewCIHBS().getBonds())

        self.outer_cihbs_obj.setParameters(self.inner_cihbs_obj.getInnerCIHBS())
        self.outer_cihbs_obj.connectResidueCIHBS()
        new_cihbs.extendBond(self.outer_cihbs_obj.getNewCIHBS().getBonds())

        mask.extra_bonds = new_cihbs.serialize()
        mask.hide_atoms = None
        mask.show_atoms = self.inner_cihbs_obj.getInnerCIHBS_serial_number()
        return mask