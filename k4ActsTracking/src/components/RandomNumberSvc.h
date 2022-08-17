// D. Elitez, August 2022


#ifndef RANDOMNUMBERSVC_H
#define RANDOMNUMBERSVC_H

#include "IRandomNumberSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/EventIDBase.h"


#include <random>

using RandomEngine = std::mt19937;

class RandomNumberSvc : public extends<Service, IRandomNumberSvc> {
public:

  RandomNumberSvc(const std::string& name, ISvcLocator* svc);

  virtual ~RandomNumberSvc();

  virtual StatusCode initialize() final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;

  virtual uint64_t getSeed();

  virtual RandomEngine getRdn();

  uint64_t generateSeed();

  RandomEngine spawnGenerator(uint64_t seed);


private:

  MsgStream m_log;

  uint64_t m_seedPtr;

  RandomEngine m_rdn;

  Gaudi::Property<int> algNum{this, "algNum", 0, "Option for propagation mode."};

  Gaudi::Property<int> evNum{this, "evNum", 0, "Option for propagation mode."};

};
inline uint64_t RandomNumberSvc::getSeed() { return m_seedPtr; }
inline RandomEngine RandomNumberSvc::getRdn() {return m_rdn; }
#endif
