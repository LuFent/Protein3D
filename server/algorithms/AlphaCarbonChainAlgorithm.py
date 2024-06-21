from server.Algorithm.Algorithm import Algorithm, Mask


class AlphaCarbonChainAlgorithm(Algorithm):
    code_name = "alpha_carbon_skeleton"
    label = "Alpha-Carbon Skeleton"
    icon = "AlphaCarbonAlgorithmLogo.svg"

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