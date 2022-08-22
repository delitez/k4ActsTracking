# D. Elitez, August 2022
# Fatras Simulation algorithm for gaudi4acts

import os
import sys
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import FatrasSimulation
from Configurables import ParticleGunAlg
from Configurables import GeoSvc

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits


algList = []

a = ParticleGunAlg("MyParticleGunAlg")
a.d0Sigma = 15 * actsUnits.um
a.z0Sigma = 55 * actsUnits.mm
a.tSigma = 1 * actsUnits.ns
a.nMultiplicity = 5;
a.nParticles = 50;
algList.append(a)

#
b = FatrasSimulation("MyFatrasSimulation")
algList.append(b)

c = GeoSvc("GeoSvc")
c.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
c.debugGeometry = True
c.outputFileName = "MyObjFileParticleGunTest"

from Configurables import THistSvc
THistSvc().Output = ["rec DATAFILE='fatrasSim_test.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().OutputLevel = DEBUG
THistSvc().PrintAll = True
THistSvc().AutoSave = True
THistSvc().AutoFlush = True

from Configurables import ApplicationMgr
ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=20, ExtSvc=[c], OutputLevel=DEBUG)
