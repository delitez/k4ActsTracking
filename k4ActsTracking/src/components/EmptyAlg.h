#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "ParticleGunAlg.h"


class EmptyAlg : public GaudiAlgorithm {
public:
  explicit EmptyAlg(const std::string&, ISvcLocator*);
  virtual ~EmptyAlg();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // member variable
  int m_member = 0;
  DataObjectHandle<AnyDataWrapper<SimParticleContainer>> p_partvec{"/Event/testVec", Gaudi::DataHandle::Reader, this};


};
