from enum import Enum

class CHIBSBond(Enum):
    physicalOperator = 1
    donorAcceptor = 2
    acceptorAcceptor = 3

class Atoms(Enum):
    singleTargetElements = ["N", "O", "S"]
    allAtomsElements = ["N", "O", "S", "C"]

    donor = "N"
    acceptor = ["O", "S"]

class Groups(Enum):
    allResidues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
                   "HIS", "ILE", "LYS", "MET", "PHE", "PRO", "SER", "THR",
                   "TPR", "TYR", "VAL"]

    sycleResidues = ["PRO", "TYR", "HIS", "TRP"]
    singleResidues = {"LYS": "NZ", "SER": "OG", "CYS": "SG", "MET": "SD", "THR": "OG1"}
    complexResidues = {"ARG": "CZ", "ASN": "CG", "GLN": "CD", "ASP": "CG", "GLU": "CD"}

    noCIHBSResidue = ["GLY", "ALA", "VAL", "LEU", "ISO"]
    residueConnectors = ["TYR", "SER", "THR"]
    antiConnectors = ["PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE"]
    doubleBondAcid = {"ASN": "OD1", "GLN": "OE1", "ARG": "NH1", "ASP": "OD1", "GLU": "OE1"}

    plusIONResidue = {"ARG": "NH1", "LYS": "NZ", "HIS": "ND1", "PRO": "N"}
    minusIONResidue = {"ASP": ["OE1", "OE2"], "GLN": ["OE1", "OE2"]}

    CIHBSInitialisators = ["PRO", "MET"]