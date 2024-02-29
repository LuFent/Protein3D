from enum import Enum

class CHIBSBond(Enum):
    physicalOperator = 1
    donorAcceptor = 2
    acceptorAcceptor = 3

class AminoMap(Enum):
    HIS = ({-1: "ND1"}, {1: "ND2"})
    TRP = ({-1: "NE1"}, {1: "NE2"})
    PRO = ({-1: "N"}, {1: "O"})
    GLY = ({-1: "N"}, {1: "O"})
    TYR = ({-1: "OH"}, {1: "OH"})

    THR = ({-1: "OG1"}, {2: "OG1"})
    SER = ({-1: "OG"}, {2: "OG"})
    CYS = ({-1: "SG1"}, {2: "SG1"})
    MET = ({-1: "SD"}, {2: "SD"})

    GLU = ({-1: "OE1"}, {1: "OE1"}, {2: "OE2"})
    ASP = ({-1: "OD1"}, {1: "OD1"}, {2: "OD2"})

    GLN = ({-2: "NE2"}, {2: "OE1"})
    ASN = ({-2: "ND2"}, {2: "OD1"})

    ARG = ({-1: "NE"}, {-2: "NH1"},{-1: "NH2"}, {1: "NH2"})

    LYS = ({-2: "NZ"}, {1: "NZ"})

    def getMap(self):
        map = [
            [None, None, None, None],
            [None,
                [AminoMap.HIS, AminoMap.TRP, AminoMap.PRO, AminoMap.GLY, AminoMap.TYR],
                [AminoMap.THR, AminoMap.SER, AminoMap.CYS, AminoMap.MET],
                [AminoMap.GLU, AminoMap.ASP]
                ],
            [None,
                None,
                [AminoMap.GLN, AminoMap.ASN],
                None,
                ],
            [None, None, None, None],
            [None,
                [AminoMap.ARG],
                None,
                None,
                None,
                ]
        ]
        return map

class Atoms(Enum):
    singleTargetElements = ["N", "O", "S"]
    allAtomsElements = ["N", "O", "S", "C"]

    donor = "N"
    acceptor = ["O", "S"]

class Groups(Enum):
    sycleResidues = ["PRO", "TYR", "HIS", "TRP"]
    singleResidues = {"LYS": "NZ", "SER": "OG", "CYS": "SG", "MET": "SD", "THR": "OG1"}
    complexResidues = {"ARG": "CZ", "ASN": "CG", "GLN": "CD", "ASP": "CG", "GLU": "CD"}

    noCIHBSResidue = ["GLY", "ALA", "VAL", "LEU", "ISO"]
    residueConnectors = ["TYR", "SER", "THR"]
    antiConnectors = ["PRO", "ALA", "VAL", "ILE", "LEU", "MET", "PHE"]
    doubleBondAcid = {"ASN": "OD1", "GLN": "OE1", "ARG": "NH1", "ASP": "OD1", "GLU": "OE1"}

    plusIONResidue = [["ARG", "LYS", "HIS", "PRO"], ["NH1", "NZ", "ND1", "N"]]
    minusIONResidue = [["ASP", "GLN"], ["OE1", "OE2"]]

    CIHBSInitialisators = ["PRO", "MET"]