import os
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import EventGeneratorAlg

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits

algList = []


a = EventGeneratorAlg("MyEventGeneratorAlg")
a.d0Sigma = 15 * actsUnits.um
a.z0Sigma = 55 * actsUnits.mm
a.tSigma = 1 * actsUnits.ns
a.nMultiplicity = 5;
a.nParticles = 10;
a.objectPath = "/Event/MyParticle1"
algList.append(a)



from Configurables import ApplicationMgr

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=1, ExtSvc=[], OutputLevel=DEBUG)
