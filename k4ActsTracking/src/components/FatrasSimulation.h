// D. Elitez, August 2022
// Fatras Simulation algorithm for gaudi4acts

#ifndef FatrasSimulation_H
#define FatrasSimulation_H

//GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/AnyDataHandle.h"
#include "GaudiKernel/AnyDataWrapper.h"
#include "GaudiKernel/DataObjectHandle.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"

//BOOST
#include <boost/container/flat_set.hpp>

//ACTS + ActsFatras
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Kernel/InteractionList.hpp"
#include "ActsFatras/Kernel/Simulation.hpp"
#include "ActsFatras/Physics/Decay/NoDecay.hpp"
#include "ActsFatras/Physics/ElectroMagnetic/PhotonConversion.hpp"
#include "ActsFatras/Physics/StandardInteractions.hpp"
#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"
#include "ActsFatras/Selectors/SurfaceSelectors.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

#include "IGeoSvc.h"
#include "TTree.h"

struct HitSurfaceSelector {
  bool sensitive = true;
  bool material  = true;
  bool passive   = false;

  bool operator()(const Acts::Surface& surface) const {
    bool isSensitive = surface.associatedDetectorElement() != nullptr;
    bool isMaterial  = surface.surfaceMaterial() != nullptr;
    bool isPassive   = not(isSensitive or isMaterial);
    return (isSensitive and sensitive) or (isMaterial and material) or (isPassive and passive);
  }
};

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

struct GeometryIdGetter {
  constexpr Acts::GeometryIdentifier operator()(Acts::GeometryIdentifier geometryId) const { return geometryId; }
  constexpr Acts::GeometryIdentifier operator()(Acts::GeometryIdentifier::Value encoded) const {
    return Acts::GeometryIdentifier(encoded);
  }
  template <typename T>
  constexpr Acts::GeometryIdentifier operator()(const std::pair<Acts::GeometryIdentifier, T>& mapItem) const {
    return mapItem.first;
  }
  template <typename T>
  inline auto operator()(const T& thing) const -> decltype(thing.geometryId(), Acts::GeometryIdentifier()) {
    return thing.geometryId();
  }
  template <typename T>
  inline auto operator()(std::reference_wrapper<T> thing) const
      -> decltype(thing.get().geometryId(), Acts::GeometryIdentifier()) {
    return thing.get().geometryId();
  }
};

struct CompareGeometryId {
  using is_transparent = void;
  template <typename Left, typename Right> constexpr bool operator()(Left&& lhs, Right&& rhs) const {
    return GeometryIdGetter()(lhs) < GeometryIdGetter()(rhs);
  }
};

using SimParticleContainer = ::boost::container::flat_set<::ActsFatras::Particle, detail::CompareParticleId>;
template <typename T> using GeometryIdMultiset = boost::container::flat_multiset<T, CompareGeometryId>;
using SimHitContainer                          = GeometryIdMultiset<::ActsFatras::Hit>;

class FatrasSimulation : public GaudiAlgorithm {
public:
  SmartIF<IGeoSvc> m_geoSvc;

  FatrasSimulation(const std::string& name, ISvcLocator* svc);

  virtual ~FatrasSimulation();

  virtual StatusCode initialize() final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;

  virtual StatusCode initializeTrees() final;

  virtual StatusCode cleanTrees() final;

private:
  ITHistSvc* m_ths{nullptr};
  TTree*     m_outputTree{nullptr};

  std::vector<float> testVar;
  std::vector<float> m_x;   ///< global x
  std::vector<float> m_y;   ///< global y
  std::vector<float> m_z;   ///< global z
  std::vector<float> m_dx;  ///< global direction x
  std::vector<float> m_dy;  ///< global direction y
  std::vector<float> m_dz;  ///< global direction z

  std::vector<int>                                       m_volumeID;     ///< volume identifier
  std::vector<int>                                       m_boundaryID;   ///< boundary identifier
  std::vector<int>                                       m_layerID;      ///< layer identifier if
  std::vector<int>                                       m_approachID;   ///< surface identifier
  std::vector<int>                                       m_sensitiveID;  ///< surface identifier
  DataObjectHandle<AnyDataWrapper<SimParticleContainer>> p_partvec{"/Event/testVec", Gaudi::DataHandle::Reader, this};
};

#endif
