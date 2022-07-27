
// D. Elitez, July 2022
// Based on eic/juggler

#ifndef PropagatorAlg_H
#define PropagatorAlg_H


#include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "IGeoSvc.h"
#include "IPropagatorAlg.h"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

#include "ProcessCode.hpp"
#include "AlgorithmContext.hpp"
#include "CommonGeometry.hpp"
#include "CommonOptions.hpp"
#include "IBaseDetector.hpp"
#include "MagneticFieldOptions.hpp"

#include "WhiteBoard.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <random>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "CommonOptions.hpp"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Options.hpp"

#include <iostream>

#include <boost/program_options.hpp>

class PropagatorInterface;

using RecordedMaterial      = Acts::MaterialInteractor::result_type;
using RecordedMaterialTrack = std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;
using PropagationOutput     = std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;
std::optional<Acts::BoundSymMatrix> generateCovariance(std::mt19937& rng, std::normal_distribution<double>& gauss);

class PropagatorAlg : public GaudiAlgorithm {
public:
  struct Config {
    /// Instance of a propagator wrapper that performs the actual propagation
    std::shared_ptr<PropagatorInterface> propagatorImpl = nullptr;

    /// how to set it up
    std::mt19937 randomNumberSvc{0};

    /// phi range
    std::pair<double, double> phiRange = {-M_PI, M_PI};
    /// eta range
    std::pair<double, double> etaRange = {-4., 4.};
    /// pt range
    std::pair<double, double> ptRange = {100 * Acts::UnitConstants::MeV, 100 * Acts::UnitConstants::GeV};

    /// The step collection to be stored
    std::string propagationStepCollection = "PropagationSteps";

    /// The material collection to be stored
    std::string propagationMaterialCollection = "RecordedMaterialTracks";

    /// The covariance values
    Acts::BoundVector covariances = Acts::BoundVector::Zero();

    /// The correlation terms
    Acts::BoundSymMatrix correlations = Acts::BoundSymMatrix::Identity();
  };

  const Config& config() const { return m_cfg; }

private:
  Config m_cfg;

  Config readPropagationConfig(const boost::program_options::variables_map& vm);

  Gaudi::Property<int> mode{this, "_mode", 0, "Option for proapgation mode."};

  Gaudi::Property<bool> sterileLogger{this, "_sterileLogger", false, "Option to switch the logger to sterile."};

  Gaudi::Property<bool> debugOutput{this, "_debugOutput", false, "Option for the debug output."};

  Gaudi::Property<bool> energyLoss{this, "_energyLoss", true, "Option to modify the behavior of the material interaction: energy loss."};

  Gaudi::Property<bool> multipleScattering{this, "_multipleScattering", true, "Option to modify the behavior of the material interaction: scattering"};

  Gaudi::Property<bool> recordMaterialInteractions{this, "_recordMaterialInteractions", true, "Option to modify the behavior of the material interaction: record"};

  Gaudi::Property<int> ntests{this, "_ntests", 0, "Option for number of particles"};

  Gaudi::Property<bool> covarianceTransport{this, "_covarianceTransport", false, "Option for covariance transport."};

  Gaudi::Property<double> d0Sigma{this, "_d0Sigma", 0, "Option for d0 gaussian sigma"};

  Gaudi::Property<double> z0Sigma{this, "_z0Sigma", 0, "Option for z0 gaussian sigma"};

  Gaudi::Property<double> phiSigma{this, "_phiSigma", 0, "Option for phi gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> thetaSigma{this, "_thetaSigma", 0, "Option for theta gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> qpSigma{this, "_qpSigma", 0, "Option for qp gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> tSigma{this, "_tSigma", 0, "Option for t gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> ptLoopers{this, "_ptLoopers", 0, "Option for looper protection"};

  Gaudi::Property<double> maxStepSize{this, "_maxStepSize", 0, "Option for Max step size steering"};

  std::random_device                  rd{};
  std::mt19937                        rng{rd()};
  std::optional<Acts::BoundSymMatrix> generateCovariance(std::mt19937& rng, std::normal_distribution<double>& gauss);

public:

  SmartIF<IGeoSvc> m_geoSvc;

  explicit PropagatorAlg(const std::string&, ISvcLocator*);

  virtual ~PropagatorAlg();

  StatusCode initialize() override final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;

};


#endif
