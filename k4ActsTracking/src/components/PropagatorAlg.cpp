// D. Elitez, July 2022

#include "PropagatorAlg.h"
#include <boost/program_options.hpp>
#include "GaudiKernel/Service.h"
#include "TGeoManager.h"
#include "TTree.h"
#include "GaudiKernel/ServiceHandle.h"
DECLARE_COMPONENT(PropagatorAlg)

using namespace Acts;

std::optional<Acts::BoundSymMatrix> PropagatorAlg::generateCovariance() {
  Rndm::Numbers gauss( randSvc(), Rndm::Gauss( 0., 1. ) );

  if (covarianceTransport) {
    Acts::BoundSymMatrix newCov(m_cfg.correlations);

    Acts::BoundVector covs_smeared = m_cfg.covariances;
    for (size_t k = 0; k < size_t(covs_smeared.size()); ++k) {
      covs_smeared[k] *= gauss();
    }

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

StatusCode PropagatorAlg::initialize() {
  m_geoSvc = service("GeoSvc");

  if (!m_geoSvc) {
    std::cout << "Unable to locate Geometry Service. " << std::endl;
    return StatusCode::FAILURE;
  }

  if (service("THistSvc", m_ths).isFailure()) {
    error() << "Couldn't get THistSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree = new TTree("hits", "PropagatorAlg hits ntuple");
  if (m_ths->regTree("/rec/NtuplesHits", m_outputTree).isFailure()) {
    error() << "Couldn't register hits tree" << endmsg;
  }
  initializeTrees();

  return StatusCode::SUCCESS;
}

StatusCode PropagatorAlg::execute() {
  cleanTrees();

  static int tSeed = 5;
  std::mt19937                     rng{tSeed++};

  // std::uniform_real_distribution<double> phiDist(0, 2 * M_PI);
  // std::uniform_real_distribution<double> etaDist(-4.0, 4.0);
  // std::uniform_real_distribution<double> ptDist(10., 20.);
  // std::uniform_real_distribution<double> qDist(0., 1.);

  Rndm::Numbers gauss( randSvc(), Rndm::Gauss( 0., 1. ) );
  Rndm::Numbers phiDist( randSvc(), Rndm::Flat( 0., 2*M_PI ) );
  Rndm::Numbers etaDist( randSvc(), Rndm::Flat( -4., 4. ) );
  Rndm::Numbers ptDist( randSvc(), Rndm::Flat( 10., 20. ) );
  Rndm::Numbers qDist( randSvc(), Rndm::Flat( 0., 1. ) );

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0., 0., 0.));

  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(ntests);

  std::vector<RecordedMaterialTrack> recordedMaterial;
  if (recordMaterialInteractions) {
    recordedMaterial.reserve(ntests);
  }
  for (size_t i = 0; i < 42; i++) {
    testVar.push_back(i);
  }

  const SimParticleContainer* p_partvect = p_partvec.get();

    for (auto i = p_partvect->begin(); i < p_partvect->end(); i++) {

    double d0     = d0Sigma * gauss();    //TODO :: set from 4pos of the particle position
    double z0     = z0Sigma * gauss();     /// parameter aus singleboundtrackparameters, see link from 10.8.
    double phi    = phiDist();           /// random numbers-randomnumbersvc?
    double eta    = etaDist();
    double theta  = 2 * atan(exp(-eta));
    double t      = tSigma * gauss();


    double pt     = i->transverseMomentum();
    double p      = pt / sin(theta);
    double charge = i->charge();
    double qop    = charge / p;


    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;


    Acts::Vector3 sPosition(0., 0., 0.);
    Acts::Vector3 sMomentum(0., 0., 0.);

    // The covariance generation
    auto cov = generateCovariance();


    auto tGeometry = m_geoSvc->trackingGeometry();

    auto                       bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2 * Acts::UnitConstants::T));
    Acts::BoundTrackParameters startParameters(surface, std::move(pars), std::move(cov));

    auto stepper     = Acts::EigenStepper<>{std::move(bField)};
    using Stepper    = std::decay_t<decltype(stepper)>;
    using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
    Acts::Navigator::Config navCfg{tGeometry};
    Acts::Navigator         navigator(navCfg);
    Propagator              propagator(std::move(stepper), std::move(navigator));

    PropagationOutput pOut;

    pOut = executeTest(propagator, startParameters);

    for (auto& step : pOut.first) {
      Acts::GeometryIdentifier::Value volumeID    = 0;
      Acts::GeometryIdentifier::Value boundaryID  = 0;
      Acts::GeometryIdentifier::Value layerID     = 0;
      Acts::GeometryIdentifier::Value approachID  = 0;
      Acts::GeometryIdentifier::Value sensitiveID = 0;

      if (step.surface) {
        auto geoID  = step.surface->geometryId();
        volumeID    = geoID.volume();
        boundaryID  = geoID.boundary();
        layerID     = geoID.layer();
        approachID  = geoID.approach();
        sensitiveID = geoID.sensitive();
      }

      if (sensitiveID >= sensitiveIDopt) {
        // a current volume overwrites the surface tagged one
        if (step.volume != nullptr) {
          volumeID = step.volume->geometryId().volume();
        }
        // now fill
        m_sensitiveID.push_back(sensitiveID);
        m_approachID.push_back(approachID);
        m_layerID.push_back(layerID);
        m_boundaryID.push_back(boundaryID);
        m_volumeID.push_back(volumeID);

        m_x.push_back(step.position.x());
        m_y.push_back(step.position.y());
        m_z.push_back(step.position.z());
        auto direction = step.momentum.normalized();
        m_dx.push_back(direction.x());
        m_dy.push_back(direction.y());
        m_dz.push_back(direction.z());
      }
    }
    m_outputTree->Fill();
  } ///// END OF FOR LOOP OVER NTESTS

  return StatusCode::SUCCESS;
}

StatusCode PropagatorAlg::finalize() { return StatusCode::SUCCESS; }

StatusCode PropagatorAlg::initializeTrees() {
  m_outputTree->Branch("testVar", &testVar);
  m_outputTree->Branch("g_x", &m_x);
  m_outputTree->Branch("g_y", &m_y);
  m_outputTree->Branch("g_z", &m_z);
  m_outputTree->Branch("d_x", &m_dx);
  m_outputTree->Branch("d_y", &m_dy);
  m_outputTree->Branch("d_z", &m_dz);
  m_outputTree->Branch("volume_id", &m_volumeID);
  m_outputTree->Branch("boundary_id", &m_boundaryID);
  m_outputTree->Branch("layer_id", &m_layerID);
  m_outputTree->Branch("approach_id", &m_approachID);
  m_outputTree->Branch("sensitive_id", &m_sensitiveID);
  return StatusCode::SUCCESS;
}

StatusCode PropagatorAlg::cleanTrees() {
  testVar.clear();
  m_x.clear();
  m_y.clear();
  m_z.clear();
  m_dx.clear();
  m_dy.clear();
  m_dz.clear();
  m_volumeID.clear();
  m_boundaryID.clear();
  m_layerID.clear();
  m_approachID.clear();
  m_sensitiveID.clear();

  return StatusCode::SUCCESS;
}
