from V_algebra import *

# Defining constants
MaxChange = 6
MaxFile = 35
StandardSoundTone = 1500
StandardSoundPause = 40
DeltaBasisNameY = 16
MinBasisMenuY = 40
MinWindowX = 200
MinWindowY = 150
DeltaIndexX = 19
DeltaIndexY = 20
MaxScreenX = 639
MaxScreenY = 479
MinX = 70
MinY = 60
MaxX = 600
MaxY = 450
DeltaFont = 8
StandardPoint = 3
BasisString = 'Basis name'
UVWString = 'Direction:'
WidthString = 'Width model [•]: '
LengthString = 'Length model[•]: '
ReflexString = 'Max reflex: '
WaitString = 'Wait '
IndexActiveString = 'Index active '
IndexPassiveString = 'Index passive '
WindowActiveString = 'Window active '
WindowPassiveString = 'Window passive '
LengthStringUVW = 11
LengthStringModel = 4
LengthStringReflex = 2
LengthStringIndex = 13
LengthStringWindow = 14
BasisX = DeltaFont * len(BasisString)
DirectionX = 20 + BasisX + DeltaFont * len(UVWString)
ModelX = 5 + DirectionX + DeltaFont * (LengthStringUVW + len(LengthString))
ReflexX = 5 + ModelX + DeltaFont * (LengthStringModel + len(ReflexString))
BasisY = 20
DirectionY = 20
ModelY = 20
ReflexY = 20


