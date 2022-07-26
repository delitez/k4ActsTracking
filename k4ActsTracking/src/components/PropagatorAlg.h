
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



//#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ProcessCode.hpp"
//#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "AlgorithmContext.hpp"
#include "CommonGeometry.hpp"
#include "CommonOptions.hpp"
#include "IBaseDetector.hpp"
#include "MagneticFieldOptions.hpp"

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
  //void write(std::ostream& os) const final;

  SmartIF<IGeoSvc> m_geoSvc;

  explicit PropagatorAlg(const std::string&, ISvcLocator*);

  virtual ~PropagatorAlg();

  StatusCode initialize() override final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;

};


#endif
