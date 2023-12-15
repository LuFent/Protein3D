from enum import Enum

class Bond(Enum):
    SINGLE = (1,"Grey")
    DOUBLE = (2, "Yellow")

class CHIBSBond(Enum):
    O = "Red"
    C = "Grey"
    N = "Green"
    R = "Blue"

class Alkaline(Enum):
    LYS = 1
    ARG = 2

class Acidic(Enum):
    ASP = 3
    GLU = 4

class Neutral(Enum):
    ASN = 5
    GLN = 6

class WeaklyPolar(Enum):
    SER = 7
    THR = 8

class NonPolar(Enum):
    GLY = 9
    VAL = 10
    LEU = 11
    Ile = 12
    ALA = 13

class SulfidCysbridges(Enum):
    CYS = 14
    Met = 15

class Cyclic(Enum):
    PRO = 16
    PHE = 17
    TYR = 18
    HIS = 19
    TRP = 20

