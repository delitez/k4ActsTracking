// D. Elitez, August 2022


#ifndef IRANDOMNUMBERSVC_H
#define IRANDOMNUMBERSVC_H

#include <GaudiKernel/IService.h>
#include <unordered_map>


class GAUDI_API IRandomNumberSvc : virtual public IService {
public:


public:
  DeclareInterfaceID(IRandomNumberSvc, 1, 0);
  virtual uint64_t getSeed() = 0;
  virtual ~IRandomNumberSvc() {}
};

#endif
