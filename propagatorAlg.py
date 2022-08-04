import os
import sys
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import PropagatorAlg
from Configurables import GeoSvc

sys.path.append('/home/delitez/ACTS/spack/k4actstracking')
import actsUnits


algList = []

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
b.sensitiveIDopt = 0

algList.append(b)


a = GeoSvc("GeoSvc")
a.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
a.debugGeometry = True
a.outputFileName = "MyObjFile"


from Configurables import ApplicationMgr

from Configurables import THistSvc
THistSvc().Output = ["rec DATAFILE='propagatorAlgOutput.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().OutputLevel = DEBUG
THistSvc().PrintAll = True
THistSvc().AutoSave = True
THistSvc().AutoFlush = True

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=500, ExtSvc=[a], OutputLevel=DEBUG)
