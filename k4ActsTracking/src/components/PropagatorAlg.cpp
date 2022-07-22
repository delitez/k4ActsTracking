// D. Elitez, July 2022

#include "PropagatorAlg.h"
#include <boost/program_options.hpp>
#include "GaudiKernel/Service.h"
#include "TGeoManager.h"

//DECLARE_COMPONENT(PropagatorAlg)

using namespace Acts;

std::optional<Acts::BoundSymMatrix> PropagatorAlg::generateCovariance(std::mt19937&                     rng,
                                                                      std::normal_distribution<double>& gauss) {
  if (m_cfg.covarianceTransport) {
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
  propagationSteps.reserve(m_cfg.ntests);

  std::vector<RecordedMaterialTrack> recordedMaterial;
  if (m_cfg.recordMaterialInteractions) {
    recordedMaterial.reserve(m_cfg.ntests);
  }

  // loop over number of particles
  for (size_t it = 0; it < m_cfg.ntests; ++it) {
    /// get the d0 and z0
    double d0     = m_cfg.d0Sigma * gauss(rng);
    double z0     = m_cfg.z0Sigma * gauss(rng);
    double phi    = phiDist(rng);
    double eta    = etaDist(rng);
    double theta  = 2 * atan(exp(-eta));
    double pt     = ptDist(rng);
    double p      = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop    = charge / p;
    double t      = m_cfg.tSigma * gauss(rng);

    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;
    // some screen output

    Acts::Vector3 sPosition(0., 0., 0.);
    Acts::Vector3 sMomentum(0., 0., 0.);

    // The covariance generation
    auto cov = generateCovariance(rng, gauss);

    // Setup and parse propagationSteps

    auto desc = ActsExamples::Options::makeDefaultOptions();
    ActsExamples::Options::addSequencerOptions(desc);
    ActsExamples::Options::addGeometryOptions(desc);
    ActsExamples::Options::addMaterialOptions(desc);
    ActsExamples::Options::addMagneticFieldOptions(desc);
    ActsExamples::Options::addRandomNumbersOptions(desc);
    //ActsExamples::Options::addPropagationOptions(desc);
    ActsExamples::Options::addOutputOptions(desc, ActsExamples::OutputFormat::Root | ActsExamples::OutputFormat::Obj);

    ActsExamples::IBaseDetector*          detector;
    boost::program_options::variables_map vm;
    auto                                  geometry          = ActsExamples::Geometry::build(vm, *detector);
    auto                                  tGeometry         = geometry.first;
    auto                                  contextDecorators = geometry.second;

    auto bField = ActsExamples::Options::readMagneticField(vm);

    bool rootOutput = vm["output-root"].template as<bool>();
    bool objOutput  = vm["output-obj"].template as<bool>();

    auto stepper     = Acts::EigenStepper<>{std::move(bField)};
    using Stepper    = std::decay_t<decltype(stepper)>;
    using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
    Acts::Navigator::Config navCfg{tGeometry};
    Acts::Navigator         navigator(navCfg);
    Propagator              propagator(std::move(stepper), std::move(navigator));

    auto randomNumberSvc = rng;
    // // Read the propagation config and create the algorithms
    auto pAlgConfig            = readPropagationConfig(vm);
    pAlgConfig.randomNumberSvc = randomNumberSvc;
    pAlgConfig.sterileLogger   = not rootOutput and not objOutput;

    // pAlgConfig.propagatorImpl =
    //     std::make_shared<PropagatorInterface::ConcretePropagator<Propagator>>(
    //         std::move(propagator));
    //
    //     sequencer.addAlgorithm(std::make_shared<ActsExamples::PropagationAlgorithm>(
    //         pAlgConfig, logLevel));
  }

  return StatusCode::SUCCESS;
}

StatusCode PropagatorAlg::execute() { return StatusCode::SUCCESS; }

StatusCode PropagatorAlg::finalize() { return StatusCode::SUCCESS; }
