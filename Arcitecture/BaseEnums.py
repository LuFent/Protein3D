from enum import Enum

class CHIBSBond(Enum):
    outer = "new"

class Atoms(Enum):
    singleTargetAtoms = ["N", "O", "S"]

class Groups(Enum):
    passive = ["ALA", "VAL", "LEU", "ILE", "PHE"]
    sycleResiduces = ["PRO", "TYR", "HIS", "TRP"]
    singleResidues = ["LYS", "CYS", "MET", "SER", "THR"]

    initiators = ["PRO", "MET"]
    valvesOneMulty = ["SER", "THR", "ASP", "GLU"]
    valvesMultyMulty = ["ASN", "GLN", "LYS", "ARG"]
    valvesOneOne = ["HIS", "TYR"]
    inverse = ["TYR", "TPR"]
