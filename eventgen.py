import os
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import ParticleGunAlg

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits

algList = []


a = ParticleGunAlg("MyParticleGunAlg")
a.d0Sigma = 15 * actsUnits.um
a.z0Sigma = 55 * actsUnits.mm
a.tSigma = 1 * actsUnits.ns
a.nMultiplicity = 5;
a.nParticles = 10;
algList.append(a)
from Configurables import EmptyAlg
d=EmptyAlg("MyEmptyAlg")
algList.append(d)


from Configurables import ApplicationMgr

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=1, ExtSvc=[], OutputLevel=DEBUG)
