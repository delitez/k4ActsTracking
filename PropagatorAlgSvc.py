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
b._mode = 0
b._sterileLogger = False
#b._debugOutput = False
b._energyLoss = True
b._multipleScattering = True
b._recordMaterialInteractions = True
b._ntests = 100
b._d0Sigma = 15 * actsUnits.um
b._z0Sigma = 55 * actsUnits.mm
b._phiSigma = 0.001
b._thetaSigma = 0.001
b._covarianceTransport = False
b._qpSigma = 0.0001 / 1 * actsUnits.GeV
b._tSigma = 1 * actsUnits.ns
b._ptLoopers = 500 * actsUnits.MeV
b._maxStepSize = 3 * actsUnits.m

algList.append(b)


a = GeoSvc("GeoSvc")
a.detectors = ["/home/delitez/ACTS/acts/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
a.debugGeometry = True
a.outputFileName = "MyObjFile"


from Configurables import ApplicationMgr

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=2, ExtSvc=[], OutputLevel=INFO)
