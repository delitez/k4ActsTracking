# Declaring ACTS constants for Gaudi Properties
import math
print("Initializing actsUnits.py")
# Length, native unit mm
fm = float(1e-12)
pm = float(1e-9)
um = float(1e-3)
nm = float(1e-6)
mm = float(1.0)
cm = float(10.0)
m = float(1e3)
km = float(1e6)
# Shortcuts for commonly used area and volume units. This intentionally
# contains not all possible combinations to avoid cluttering the namespace.
# Missing area or volume units can always be defined on the fly using the
# existing length units e.g. 1fm³ -> 1fm * 1fm * 1fm
# Area, native unit mm²
mm2 = mm * mm
cm2 = cm * cm
m2 = m * m
# Volume, native unit mm³
mm3 = mm * mm * mm
cm3 = cm * cm * cm
m3 = m * m * m
# Time, native unit mm = [speed-of-light * time] = mm/s * s
s = float(299792458000.0)
fs = float(1e-15) * s
ps = float(1e-12) * s
ns = float(1e-9) * s
us = float(1e-6) * s
ms = float(1e-3) * s
min = float(60.0) * s
h = float(3600.0) * s
# Angles, native unit radian
mrad = float(1e-3)
rad = float(1.0)
degree = float(0.017453292519943295)  # pi / 180
# Energy/mass/momentum, native unit GeV
eV = float(1e-9)
keV = float(1e-6)
MeV = float(1e-3)
GeV = float(1.0)
TeV = float(1e3)
J = float(6241509074.460763) * GeV
# atomic mass unit u
u = float(0.93149410242)
#     1eV/c² == 1.782662e-36kg
#    1GeV/c² == 1.782662e-27kg
# ->     1kg == (1/1.782662e-27)GeV/c²
# ->      1g == (1/(1e3*1.782662e-27))GeV/c²
g = float(1.0) / float(1.782662e-24)
kg = float(1.0) / float(1.782662e-27)
# Charge, native unit e (elementary charge)
e = float(1.0)
C = J / eV
# Magnetic field, native unit GeV/(e*mm)
T = float(0.000299792458)  # equivalent to c in appropriate SI units
Gauss = float(1e-4) * T
kGauss = float(1e-1) * T
# Amount of substance, native unit mol
mol = float(1.0)
print("Acts units are initialized from actsUnits.py")
