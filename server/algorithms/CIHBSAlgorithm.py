from server.Algorithm.Arcitecture import CIHBS, CIHBS_selector
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
        self.cihbs_selector_obj = CIHBS_selector.SelectorCIHBS()

    def _execute(self, structure) -> Mask:
        self.set_cihbs()
        mask = Mask()
        self.cihbs_obj.setParameters(structure=structure)
        new_cihbs = self.cihbs_obj.getNewCIHBS()
        print("new_cihbs", new_cihbs)
        self.inner_cihbs_obj.checkInnerGroups(structure, self.cihbs_obj.getMainChain())

        self.physical_cihbs_obj.setParameters(args=(self.inner_cihbs_obj.getInnerCIHBSDict(),
                                                    self.inner_cihbs_obj.getInnerCIHBS(),
                                                    self.cihbs_obj.getCarbonChain()))
        self.physical_cihbs_obj.connectPhysicalOperators()
        new_cihbs.extend(self.physical_cihbs_obj.getNewCIHBS())

        self.outer_cihbs_obj.setParameters(self.inner_cihbs_obj.getInnerCIHBS())
        self.outer_cihbs_obj.connectResidueCIHBS()
        new_cihbs.extend(self.outer_cihbs_obj.getNewCIHBS())

        #!!
        self.cihbs_selector_obj.setParameters(new_cihbs)
        mask.extra_bonds = self.cihbs_selector_obj.get_selector_map()
        mask.hide_atoms = None
        mask.show_atoms = self.inner_cihbs_obj.getInnerCIHBS_serial_number()
        return mask