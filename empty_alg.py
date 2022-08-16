import os
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import EmptyAlg
from Configurables import objectTestFind
from Configurables import objectTest
from Configurables import ParticleGunAlg
from Configurables import GeoSvc
from Configurables import RandomNumberSvc
from Configurables import EventCounter
#
# algList = []
#
#
# a = EmptyAlg("MyEmptyAlg")
#
# b = objectTest("MyObjectTest")
# #algList.append(b)
# c = objectTestFind("MyObjectTestFind")
# #algList.append(c)
# d = ParticleGunAlg("MyEvGenAlg")
# algList.append(d)
# #algList.append(a)

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits

algList = []


a = ParticleGunAlg("MyParticleGunAlg")
a.d0Sigma = 15 * actsUnits.um
a.z0Sigma = 55 * actsUnits.mm
a.tSigma = 1 * actsUnits.ns
a.nMultiplicity = 5;
a.nParticles = 10;
#a.objectPath = "/Event/MyParticle4"
#algList.append(a)

b = RandomNumberSvc("MyRndNbrSvc")




d=EmptyAlg("MyEmptyAlg")
#algList.append(d)

c = GeoSvc("GeoSvc")
c.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
c.debugGeometry = True
c.outputFileName = "MyObjFileTESTparticleGun"

from Configurables import ApplicationMgr
#print("CRASH TEST")
ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=1, ExtSvc=[b], OutputLevel=DEBUG)
