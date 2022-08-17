
// D. Elitez, July 2022
// Based on eic/juggler

#ifndef PropagatorAlg_H
#define PropagatorAlg_H

#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "IGeoSvc.h"
#include "IPropagatorAlg.h"
#include "ParticleGunAlg.h"
#include<boost/container/flat_set.hpp>


#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
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
#include "Acts/Utilities/Logger.hpp"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"


#include "AlgorithmContext.hpp"
#include "CommonGeometry.hpp"
#include "CommonOptions.hpp"
#include "IBaseDetector.hpp"
#include "MagneticFieldOptions.hpp"
#include "ProcessCode.hpp"

#include "WhiteBoard.hpp"

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <random>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "CommonOptions.hpp"
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ServiceHandle.h"
#include "Options.hpp"

#include "GaudiKernel/ITHistSvc.h"

#include "TGraph.h"
#include "TH1F.h"

#include <iostream>

#include <boost/program_options.hpp>

using RecordedMaterial      = Acts::MaterialInteractor::result_type;
using RecordedMaterialTrack = std::pair<std::pair<Acts::Vector3, Acts::Vector3>, RecordedMaterial>;
using PropagationOutput     = std::pair<std::vector<Acts::detail::Step>, RecordedMaterial>;
std::optional<Acts::BoundSymMatrix> generateCovariance();

class PropagatorAlg : public GaudiAlgorithm {
public:
  struct Config {
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
  Config                              m_cfg;

  std::vector<PropagationOutput>      testvec;
  std::mt19937                        rng;
  std::optional<Acts::BoundSymMatrix> generateCovariance();



  ITHistSvc* m_ths{nullptr};
  TTree*     m_outputTree{nullptr};

  std::vector<float> testVar;
  std::vector<float> m_x;   ///< global x
  std::vector<float> m_y;   ///< global y
  std::vector<float> m_z;   ///< global z
  std::vector<float> m_dx;  ///< global direction x
  std::vector<float> m_dy;  ///< global direction y
  std::vector<float> m_dz;  ///< global direction z

  std::vector<int> m_volumeID;     ///< volume identifier
  std::vector<int> m_boundaryID;   ///< boundary identifier
  std::vector<int> m_layerID;      ///< layer identifier if
  std::vector<int> m_approachID;   ///< surface identifier
  std::vector<int> m_sensitiveID;  ///< surface identifier

  Config readPropagationConfig(const boost::program_options::variables_map& vm);

  DataObjectHandle<AnyDataWrapper<SimParticleContainer>> p_partvec{"/Event/testVec", Gaudi::DataHandle::Reader, this};

  Gaudi::Property<int> mode{this, "mode", 0, "Option for propagation mode."};

  Gaudi::Property<bool> sterileLogger{this, "sterileLogger", false, "Option to switch the logger to sterile."};

  Gaudi::Property<bool> debugOutput{this, "debugOutput", false, "Option for the debug output."};

  Gaudi::Property<bool> energyLoss{this, "energyLoss", true,
                                   "Option to modify the behavior of the material interaction: energy loss."};

  Gaudi::Property<bool> multipleScattering{this, "multipleScattering", true,
                                           "Option to modify the behavior of the material interaction: scattering"};

  Gaudi::Property<bool> recordMaterialInteractions{this, "recordMaterialInteractions", true,
                                                   "Option to modify the behavior of the material interaction: record"};

  Gaudi::Property<int> ntests{this, "ntests", 0, "Option for number of particles"};

  Gaudi::Property<bool> covarianceTransport{this, "covarianceTransport", false, "Option for covariance transport."};

  Gaudi::Property<double> d0Sigma{this, "d0Sigma", 0, "Option for d0 gaussian sigma"};

  Gaudi::Property<double> z0Sigma{this, "z0Sigma", 0, "Option for z0 gaussian sigma"};

  Gaudi::Property<double> phiSigma{this, "phiSigma", 0,
                                   "Option for phi gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> thetaSigma{this, "thetaSigma", 0,
                                     "Option for theta gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> qpSigma{this, "qpSigma", 0, "Option for qp gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> tSigma{this, "tSigma", 0, "Option for t gaussian sigma (used for covariance transport)"};

  Gaudi::Property<double> ptLoopers{this, "ptLoopers", 0, "Option for looper protection"};

  Gaudi::Property<double> maxStepSize{this, "maxStepSize", 0, "Option for Max step size steering"};

  Gaudi::Property<int> sensitiveIDopt{this, "sensitiveIDopt", 0, "Option for sensitiveID"};

  template <typename parameters_t>
  PropagationOutput executeTest(Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>& propagator,
                                parameters_t&                                            startParameters) {
    PropagationOutput pOutput;
    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Propagation Logger", Acts::Logging::INFO));

    if (mode == 0) {
      using MaterialInteractor = Acts::MaterialInteractor;
      using SteppingLogger     = Acts::detail::SteppingLogger;
      using EndOfWorld         = Acts::EndOfWorldReached;

      using ActionList        = Acts::ActionList<SteppingLogger, MaterialInteractor>;
      using AbortList         = Acts::AbortList<EndOfWorld>;
      using PropagatorOptions = Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

      const Acts::GeometryContext      geoContext;
      const Acts::MagneticFieldContext magFieldContext;

      PropagatorOptions options(geoContext, magFieldContext, Acts::LoggerWrapper{logger()});

      options.pathLimit = std::numeric_limits<double>::max();

      options.loopProtection = (startParameters.transverseMomentum() < ptLoopers);

      auto& mInteractor              = options.actionList.get<MaterialInteractor>();
      mInteractor.multipleScattering = multipleScattering;
      mInteractor.energyLoss         = energyLoss;
      mInteractor.recordInteractions = recordMaterialInteractions;

      auto& sLogger   = options.actionList.get<SteppingLogger>();
      sLogger.sterile = sterileLogger;

      options.maxStepSize = maxStepSize;

      auto result = propagator.propagate(startParameters, options);
      if (result.ok()) {
        const auto& resultValue     = result.value();
        auto        steppingResults = resultValue.template get<SteppingLogger::result_type>();

        pOutput.first = std::move(steppingResults.steps);

        if (recordMaterialInteractions) {
          auto materialResult = resultValue.template get<MaterialInteractor::result_type>();
          pOutput.second      = std::move(materialResult);
        }
      }
      return pOutput;
    }
  }

public:
  SmartIF<IGeoSvc> m_geoSvc;

//  SmartIF<IRandomNumberSvc> m_rndSvc;

  explicit PropagatorAlg(const std::string&, ISvcLocator*);

  virtual ~PropagatorAlg();

  StatusCode initialize() override final;

  virtual StatusCode execute() final;

  virtual StatusCode finalize() final;

  virtual StatusCode initializeTrees() final;

  virtual StatusCode cleanTrees() final;

};  // class

#endif
