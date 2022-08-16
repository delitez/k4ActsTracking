import os
from pprint import pprint
from Gaudi.Configuration import *

from Configurables import objectTestFind

algList = []


a = objectTestFind("MyobjectTestFind")
algList.append(a)

from Configurables import ApplicationMgr

ApplicationMgr(TopAlg=algList, EvtSel="NONE", EvtMax=1, ExtSvc=[], OutputLevel=DEBUG)
