
// D. Elitez, July 2022
// Based on eic/juggler

#ifndef PropagatorAlg_H
#define PropagatorAlg_H

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

//#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ProcessCode.hpp"
//#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "AlgorithmContext.hpp"
#include "CommonGeometry.hpp"
#include "CommonOptions.hpp"
#include "IBaseDetector.hpp"
#include "MagneticFieldOptions.hpp"
#include "PropagationOptions.hpp"
#include "WhiteBoard.hpp"

//#include "PropagationOptions.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <random>

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
//#include "IPropagatorAlg.h"
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

    /// proapgation mode
    int mode = 0;
    /// Switch the logger to sterile
    bool sterileLogger = false;
    /// debug output
    bool debugOutput = false;
    /// Modify the behavior of the material interaction: energy loss
    bool energyLoss = true;
    /// Modify the behavior of the material interaction: scattering
    bool multipleScattering = true;
    /// Modify the behavior of the material interaction: record
    bool recordMaterialInteractions = true;

    /// number of particles
    size_t ntests = 100;
    /// d0 gaussian sigma
    double d0Sigma = 15 * Acts::UnitConstants::um;
    /// z0 gaussian sigma
    double z0Sigma = 55 * Acts::UnitConstants::mm;
    /// phi gaussian sigma (used for covariance transport)
    double phiSigma = 0.001;
    /// theta gaussian sigma (used for covariance transport)
    double thetaSigma = 0.001;
    /// qp gaussian sigma (used for covariance transport)
    double qpSigma = 0.0001 / 1 * Acts::UnitConstants::GeV;
    /// t gaussian sigma (used for covariance transport)
    double tSigma = 1 * Acts::UnitConstants::ns;
    /// phi range
    std::pair<double, double> phiRange = {-M_PI, M_PI};
    /// eta range
    std::pair<double, double> etaRange = {-4., 4.};
    /// pt range
    std::pair<double, double> ptRange = {100 * Acts::UnitConstants::MeV, 100 * Acts::UnitConstants::GeV};
    /// looper protection
    double ptLoopers = 500 * Acts::UnitConstants::MeV;

    /// Max step size steering
    double maxStepSize = 3 * Acts::UnitConstants::m;

    /// The step collection to be stored
    std::string propagationStepCollection = "PropagationSteps";

    /// The material collection to be stored
    std::string propagationMaterialCollection = "RecordedMaterialTracks";

    /// covariance transport
    bool covarianceTransport = false;

    /// The covariance values
    Acts::BoundVector covariances = Acts::BoundVector::Zero();

    /// The correlation terms
    Acts::BoundSymMatrix correlations = Acts::BoundSymMatrix::Identity();
  };

  const Config& config() const { return m_cfg; }

private:
  Config m_cfg;
  Config readPropagationConfig(const boost::program_options::variables_map& vm);

  //MsgStream m_log;

  std::random_device                  rd{};
  std::mt19937                        rng{rd()};
  std::optional<Acts::BoundSymMatrix> generateCovariance(std::mt19937& rng, std::normal_distribution<double>& gauss);

public:
  explicit PropagatorAlg(const std::string&, ISvcLocator*);

  virtual ~PropagatorAlg();

  StatusCode initialize() override final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;
};

inline PropagatorAlg::Config readPropagationConfig(const boost::program_options::variables_map& vm) {
  using namespace Acts::UnitLiterals;
  PropagatorAlg::Config pAlgConfig;

  auto iphir = vm["prop-phi-range"].template as<ActsExamples::Options::Reals<2>>();
  auto ietar = vm["prop-eta-range"].template as<ActsExamples::Options::Reals<2>>();
  auto iptr  = vm["prop-pt-range"].template as<ActsExamples::Options::Reals<2>>();

  /// Material interaction behavior
  pAlgConfig.energyLoss                 = vm["prop-energyloss"].template as<bool>();
  pAlgConfig.multipleScattering         = vm["prop-scattering"].template as<bool>();
  pAlgConfig.recordMaterialInteractions = vm["prop-record-material"].template as<bool>();

  /// Create the config for the Extrapoaltion algorithm
  pAlgConfig.debugOutput = vm["prop-debug"].template as<bool>();
  pAlgConfig.ntests      = vm["prop-ntests"].template as<size_t>();
  pAlgConfig.mode        = vm["prop-mode"].template as<int>();
  pAlgConfig.d0Sigma     = vm["prop-d0-sigma"].template as<double>() * 1_mm;
  pAlgConfig.z0Sigma     = vm["prop-z0-sigma"].template as<double>() * 1_mm;
  pAlgConfig.phiSigma    = vm["prop-phi-sigma"].template as<double>();
  pAlgConfig.thetaSigma  = vm["prop-theta-sigma"].template as<double>();
  pAlgConfig.qpSigma     = vm["prop-qp-sigma"].template as<double>() / 1_GeV;
  pAlgConfig.tSigma      = vm["prop-t-sigma"].template as<double>() * 1_ns;

  pAlgConfig.phiRange    = {iphir[0], iphir[1]};
  pAlgConfig.etaRange    = {ietar[0], ietar[1]};
  pAlgConfig.ptRange     = {iptr[0] * 1_GeV, iptr[1] * 1_GeV};
  pAlgConfig.ptLoopers   = vm["prop-pt-loopers"].template as<double>() * 1_GeV;
  pAlgConfig.maxStepSize = vm["prop-max-stepsize"].template as<double>() * 1_mm;

  pAlgConfig.propagationStepCollection     = vm["prop-step-collection"].template as<std::string>();
  pAlgConfig.propagationMaterialCollection = vm["prop-material-collection"].template as<std::string>();

  /// The covariance transport
  if (vm["prop-cov"].template as<bool>()) {
    /// Set the covariance transport to true
    pAlgConfig.covarianceTransport = true;
    /// Set the covariance matrix
    pAlgConfig.covariances(Acts::BoundIndices::eBoundLoc0)   = pAlgConfig.d0Sigma * pAlgConfig.d0Sigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundLoc1)   = pAlgConfig.z0Sigma * pAlgConfig.z0Sigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundPhi)    = pAlgConfig.phiSigma * pAlgConfig.phiSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundTheta)  = pAlgConfig.thetaSigma * pAlgConfig.thetaSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundQOverP) = pAlgConfig.qpSigma * pAlgConfig.qpSigma;
    pAlgConfig.covariances(Acts::BoundIndices::eBoundTime)   = pAlgConfig.tSigma * pAlgConfig.tSigma;

    // Only if they are properly defined, assign off-diagonals
    if (vm.count("prop-corr-offd")) {
      auto readOffd = vm["prop-corr-offd"].template as<ActsExamples::Options::Reals<15>>();
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc1)    = readOffd[0];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundPhi)     = readOffd[1];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundTheta)   = readOffd[2];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundQOverP)  = readOffd[3];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundTime)    = readOffd[4];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundPhi)     = readOffd[5];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundTheta)   = readOffd[6];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundQOverP)  = readOffd[7];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundTime)    = readOffd[8];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundTheta)    = readOffd[9];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundQOverP)   = readOffd[10];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundTime)     = readOffd[11];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundQOverP) = readOffd[12];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundTheta, Acts::BoundIndices::eBoundTime)   = readOffd[13];
      pAlgConfig.correlations(Acts::BoundIndices::eBoundQOverP, Acts::BoundIndices::eBoundTime)  = readOffd[14];
    } else {
      /// Some pre-defined values (non-trivial helical correlations)
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundPhi)    = -0.8;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundQOverP) = -0.3;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundTheta)  = -0.8;
      pAlgConfig.correlations(Acts::BoundIndices::eBoundPhi, Acts::BoundIndices::eBoundQOverP)  = 0.4;
    }
  }

  return pAlgConfig;
}

#endif
