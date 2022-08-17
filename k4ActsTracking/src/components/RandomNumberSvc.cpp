// D. Elitez, August 2022


#include "RandomNumberSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "IRandomNumberSvc.h"

using namespace Gaudi;

DECLARE_COMPONENT(RandomNumberSvc)

RandomNumberSvc::RandomNumberSvc(const std::string& name, ISvcLocator* svc) : base_class(name, svc), m_log(msgSvc(), name) {}

RandomNumberSvc::~RandomNumberSvc() {}

StatusCode RandomNumberSvc::initialize() {

  auto m_seedPtr = generateSeed();

  auto m_rdn = spawnGenerator(generateSeed());

  return StatusCode::SUCCESS;
}

StatusCode RandomNumberSvc::execute() {

  SmartDataPtr<Event> evt(eventSvc(), "/Event");

  return StatusCode::SUCCESS;
}

StatusCode RandomNumberSvc::finalize() {

  return StatusCode::SUCCESS;
}

uint64_t RandomNumberSvc::generateSeed() {
uint64_t seed = 1234567890u;
const uint64_t k1 = algNum;
const uint64_t k2 = evNum;

const uint64_t id = (k1+k2) * (k1+k2+1) / 2 + k2;
return seed + id;

}

RandomEngine spawnGenerator(uint64_t seed) {
  return RandomEngine(seed);
}
