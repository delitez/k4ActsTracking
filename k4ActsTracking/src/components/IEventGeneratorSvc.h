// D. Elitez, July 2022
// Based on eic/juggler

#ifndef IEVENTGENERATORSVC_H
#define IEVENTGENERATORSVC_H

#include <GaudiKernel/IService.h>



class GAUDI_API IEventGeneratorSvc : virtual public IService {
public:


public:
  DeclareInterfaceID(IEventGeneratorSvc, 1, 0);

  virtual ~IEventGeneratorSvc() {}
};

#endif  // IGEOSVC_H
