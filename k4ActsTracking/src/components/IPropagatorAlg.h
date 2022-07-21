// D. Elitez, July 2022
// Based on eic/juggler

#ifndef IPropagatorAlg_H
#define IPropagatorAlg_H

#include <GaudiKernel/IService.h>
#include <unordered_map>



class GAUDI_API IPropagatorAlg : virtual public IService {
public:


public:
  DeclareInterfaceID(IPropagatorAlg, 1, 0);



  virtual ~IPropagatorAlg() {}
};

#endif  // IPropagatorAlg_H
