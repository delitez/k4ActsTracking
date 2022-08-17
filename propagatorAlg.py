import os
import sys
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import objectTest
from Configurables import PropagatorAlg
from Configurables import GeoSvc
from Configurables import ParticleGunAlg
from Configurables import RandomNumberSvc
from Configurables import EventCounter

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits


algList = []

# a = objectTest("MyobjectTest")
# algList.append(a)


a = ParticleGunAlg("MyParticleGunAlg")
a.d0Sigma = 15 * actsUnits.um
a.z0Sigma = 55 * actsUnits.mm
a.tSigma = 1 * actsUnits.ns
a.nMultiplicity = 5;
a.nParticles = 10;
#a.objectPath = "/Event/MyParticle4"
algList.append(a)

b = PropagatorAlg("PropagatorAlg")
#
b.mode = 0
b.sterileLogger = False
b.debugOutput = False
b.energyLoss = True
b.multipleScattering = True
b.recordMaterialInteractions = True
b.ntests = 100
b.d0Sigma = 15 * actsUnits.um
b.z0Sigma = 55 * actsUnits.mm
b.phiSigma = 0.001
b.thetaSigma = 0.001
b.covarianceTransport = False
b.qpSigma = 0.0001 / 1 * actsUnits.GeV
b.tSigma = 1 * actsUnits.ns
b.ptLoopers = 500 * actsUnits.MeV
b.maxStepSize = 3 * actsUnits.m
b.sensitiveIDopt = 1

algList.append(b)



c = GeoSvc("GeoSvc")
c.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
c.debugGeometry = True
c.outputFileName = "MyObjFileParticleGunTest"

d = RandomNumberSvc("MyRndNbrSvc")
d.algNum = len(algList);
d.evNum = b.ntests;


e = EventCounter("MyEventCounter")
algList.append(e)


from Configurables import ApplicationMgr

from Configurables import THistSvc
THistSvc().Output = ["rec DATAFILE='propagatorAlgOutput_PARTICLEGUN.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().OutputLevel = DEBUG
THistSvc().PrintAll = True
THistSvc().AutoSave = True
THistSvc().AutoFlush = True

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=20, ExtSvc=[c,d], OutputLevel=DEBUG)
