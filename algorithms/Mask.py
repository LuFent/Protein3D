class Mask:
    hide_atoms_field_name = "hideAtoms"
    show_atoms_field_name = "showAtoms"
    extra_bonds_field_name = "extraBonds"

    def __init__(self, hide_atoms=list(), show_atoms=list(), extra_bonds=dict()):
        self.hide_atoms = hide_atoms
        self.show_atoms = show_atoms
        self.extra_bonds = extra_bonds

    def flush(self):
        self.hide_atoms = list()
        self.show_atoms = list()
        self.extra_bonds = dict()

    def serialize(self):
        return {
            self.hide_atoms_field_name: self.hide_atoms,
            self.show_atoms_field_name: self.show_atoms,
            self.extra_bonds_field_name: self.extra_bonds
        }
