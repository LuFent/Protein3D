from enum import Enum

class CHIBSBond(Enum):
    physicalOperator = 1
    donorAcceptor = 2
    acceptorAcceptor = 3

class AcidGroups(Enum):
    alkaline = {1: ["LYS", "ARG"]}
    acidic = {2: ["ASP", "GLU"]}
    neutral = {3: ["ASN", "GLN"]}
    weaklyPolar = {4: ["SER", "THR"]}
    nonPolar = {5: ["GLY", "VAL", "LEU", "ILE", "ALA"]}
    sulfid = {6: ["CYS", "MET"]}
    cyclic = {7: ["PRO", "PHE", "TYR", "HIS", "TRP"]}

class Atoms(Enum):
    singleTargetElements = ["N", "O", "S"]
    allAtomsElements = ["N", "O", "S", "C"]

    donor = "N"
    acceptor = ["O", "S"]

class Groups(Enum):
    sycleResidues = ["PRO", "TYR", "HIS", "TRP"]
    singleResidues = {"LYS": "NZ", "SER": "OG", "CYS": "SG", "MET": "SD", "THR": "OG1"}
    complexResidues = {"ARG": "CZ", "ASN": "CG", "GLN": "CD", "ASP": "CG", "GLU": "CD"}

    antiConnectors = ["PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE"]
    doubleBondAcid = {"ASN": "OD1", "GLN": "OE1", "ARG": "NH1", "ASP": "OD1", "GLU": "OE1"}

    plusIONResidue = {"ARG": "NH1", "LYS": "NZ", "HIS": "ND1", "PRO": "N"}
    minusIONResidue = {"ASP": ["OE1", "OE2"], "GLN": ["OE1", "OE2"]}