// D. Elitez, August 2022
// Particle gun for gaudi4acts

#ifndef ParticleGunAlg_H
#define ParticleGunAlg_H

//GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/AnyDataHandle.h"
#include "GaudiKernel/AnyDataWrapper.h"
#include "GaudiKernel/DataObjectHandle.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

//ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

//BOOST
#include <boost/container/flat_set.hpp>

#include <cmath>
#include <random>

namespace detail {
  struct CompareParticleId {
    using is_transparent = void;
    constexpr bool operator()(const ActsFatras::Particle& lhs, const ActsFatras::Particle& rhs) const {
      return lhs.particleId() < rhs.particleId();
    }
    constexpr bool operator()(ActsFatras::Barcode lhs, const ActsFatras::Particle& rhs) const {
      return lhs < rhs.particleId();
    }
    constexpr bool operator()(const ActsFatras::Particle& lhs, ActsFatras::Barcode rhs) const {
      return lhs.particleId() < rhs;
    }
  };
}  // namespace detail

using SimParticleContainer = ::boost::container::flat_set<::ActsFatras::Particle, detail::CompareParticleId>;

class ParticleGunAlg : public GaudiAlgorithm {
public:
  SimParticleContainer genVertexParticles();

  SimParticleContainer particles;

private:
  DataObjectHandle<AnyDataWrapper<SimParticleContainer>> m_partvec{"/Event/testVec", Gaudi::DataHandle::Writer, this};

public:
  Gaudi::Property<double> d0Sigma{this, "d0Sigma", 0, "Option for d0 gaussian sigma"};

  Gaudi::Property<double> z0Sigma{this, "z0Sigma", 0, "Option for z0 gaussian sigma"};

  Gaudi::Property<double> phiSigma{this, "phiSigma", 0,
                                   "Option for phi gaussian sigma (used for covariance transport)"};
  Gaudi::Property<double> tSigma{this, "tSigma", 0, "Option for t gaussian sigma (used for covariance transport)"};

  Gaudi::Property<int> nParticles{this, "nParticles", 0, "Number of particles."};

  Gaudi::Property<int> nMultiplicity{this, "nMultiplicity", 0, "Number of multiplicity."};

  Gaudi::Property<std::string> objectPath{this, "objectPath", " ", "Path for the object."};

  ParticleGunAlg(const std::string& name, ISvcLocator* svc);

  virtual ~ParticleGunAlg();

  virtual StatusCode initialize() final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;
};

#endif
