// D. Elitez, July 2022

#include "PropagatorAlg.h"
#include <boost/program_options.hpp>
#include "GaudiKernel/Service.h"
#include "TGeoManager.h"
//#include "PropagatorInterface.hpp"
//#include "PropagationOptions.hpp"



DECLARE_COMPONENT(PropagatorAlg)

using namespace Acts;

std::optional<Acts::BoundSymMatrix> PropagatorAlg::generateCovariance(std::mt19937&                     rng,
                                                                      std::normal_distribution<double>& gauss) {
  if (covarianceTransport) {
    // We start from the correlation matrix
    Acts::BoundSymMatrix newCov(m_cfg.correlations);
    // Then we draw errors according to the error values
    Acts::BoundVector covs_smeared = m_cfg.covariances;
    for (size_t k = 0; k < size_t(covs_smeared.size()); ++k) {
      covs_smeared[k] *= gauss(rng);
    }
    // and apply a double loop
    for (size_t i = 0; i < size_t(newCov.rows()); ++i) {
      for (size_t j = 0; j < size_t(newCov.cols()); ++j) {
        (newCov)(i, j) *= covs_smeared[i];
        (newCov)(i, j) *= covs_smeared[j];
      }
    }
    return newCov;
  }
  return std::nullopt;
}

using namespace Gaudi;

PropagatorAlg::PropagatorAlg(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

PropagatorAlg::~PropagatorAlg() {}

//StatusCode PropagatorAlg::initialize(ActsExamples::AlgorithmContext& context) {
StatusCode PropagatorAlg::initialize() {
  //  m_log << MSG::INFO << "this is propagator service"<< endmsg;

  m_geoSvc = service("GeoSvc");


  if (!m_geoSvc) {
  std::cout << "Unable to locate Geometry Service. " << std::endl;
  return StatusCode::FAILURE;

}


  std::mt19937                     rng{1};
  std::normal_distribution<double> gauss(0., 1.);
  std::round(gauss(rng));

  // Setup random number distributions for some quantities
  std::uniform_real_distribution<double> phiDist(3., 5.);    // m_cfg.phiRange.first, m_cfg.phiRange.second
  std::uniform_real_distribution<double> etaDist(1.0, 2.0);  //m_cfg.etaRange.first,m_cfg.etaRange.second
  std::uniform_real_distribution<double> ptDist(10., 20.);   // m_cfg.ptRange.first,m_cfg.ptRange.second
  std::uniform_real_distribution<double> qDist(0., 1.);

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0., 0., 0.));

  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(ntests);

  std::vector<RecordedMaterialTrack> recordedMaterial;
  if (recordMaterialInteractions) {
    recordedMaterial.reserve(ntests);
  }

  // loop over number of particles
  for (size_t it = 0; it < ntests; ++it) {
    /// get the d0 and z0
    double d0     = d0Sigma * gauss(rng);
    double z0     = z0Sigma * gauss(rng);
    double phi    = phiDist(rng);
    double eta    = etaDist(rng);
    double theta  = 2 * atan(exp(-eta));
    double pt     = ptDist(rng);
    double p      = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop    = charge / p;
    double t      = tSigma * gauss(rng);

    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;
    // some screen output

    Acts::Vector3 sPosition(0., 0., 0.);
    Acts::Vector3 sMomentum(0., 0., 0.);

    // The covariance generation
    auto cov = generateCovariance(rng, gauss);

    // Setup and parse propagationSteps
    ActsExamples::IBaseDetector*          detector;

    auto tGeometry = m_geoSvc->trackingGeometry();

    auto bField = std::make_shared<ConstantBField>(Vector3(0.,0.,2 * Acts::UnitConstants::T));
    Acts::BoundTrackParameters startParameters(surface, std::move(pars),
                                                 std::move(cov));

    auto stepper     = Acts::EigenStepper<>{std::move(bField)};
    using Stepper    = std::decay_t<decltype(stepper)>;
    using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
    Acts::Navigator::Config navCfg{tGeometry};
    Acts::Navigator         navigator(navCfg);
    Propagator              propagator(std::move(stepper), std::move(navigator));

if ( mode == 0 ) {

    PropagationOutput pOutput;
    ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Propagation Logger", Acts::Logging::INFO));


    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
    using AbortList = Acts::AbortList<EndOfWorld>;
    using PropagatorOptions = Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

    const Acts::GeometryContext geoContext;
    const Acts::MagneticFieldContext magFieldContext;

    PropagatorOptions options(geoContext, magFieldContext, Acts::LoggerWrapper{logger()});

    options.pathLimit = std::numeric_limits<double>::max();

    options.loopProtection = (startParameters.transverseMomentum() < ptLoopers);

                        // Switch the material interaction on/off & eventually into logging mode
    auto& mInteractor = options.actionList.get<MaterialInteractor>();
    mInteractor.multipleScattering = multipleScattering;
    mInteractor.energyLoss = energyLoss;
    mInteractor.recordInteractions = recordMaterialInteractions;

                        // Switch the logger to sterile, e.g. for timing checks
    auto& sLogger = options.actionList.get<SteppingLogger>();
    sLogger.sterile = sterileLogger;
                        // Set a maximum step size
    options.maxStepSize = maxStepSize;

    auto result = propagator.propagate(startParameters, options);
    if (result.ok()) {
  const auto& resultValue = result.value();
  auto steppingResults =
      resultValue.template get<SteppingLogger::result_type>();

  // Set the stepping result
  pOutput.first = std::move(steppingResults.steps);
  // Also set the material recording result - if configured
  if (recordMaterialInteractions) {
    auto materialResult =
        resultValue.template get<MaterialInteractor::result_type>();
    pOutput.second = std::move(materialResult);
  }

}


   //pOutput return function

  }

}

  return StatusCode::SUCCESS;
}

StatusCode PropagatorAlg::execute() {
return StatusCode::SUCCESS;

}

StatusCode PropagatorAlg::finalize() { return StatusCode::SUCCESS; }
