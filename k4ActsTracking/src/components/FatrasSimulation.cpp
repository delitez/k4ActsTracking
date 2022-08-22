// D. Elitez, August 2022
// Fatras simulation algorithm for gaudi4acts

#include "FatrasSimulation.h"

DECLARE_COMPONENT(FatrasSimulation)

FatrasSimulation::FatrasSimulation(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

FatrasSimulation::~FatrasSimulation() {}

StatusCode FatrasSimulation::initialize() {
  MsgStream log(msgSvc(), name());

  m_geoSvc = service("GeoSvc");

  if (!m_geoSvc) {
    log << MSG::ERROR << "Unable to locate Geometry Service." << endmsg;
    return StatusCode::FAILURE;
  }

  if (service("THistSvc", m_ths).isFailure()) {
    log << MSG::ERROR << "Couldn't get THistSvc." << endmsg;
    return StatusCode::FAILURE;
  }

  m_outputTree = new TTree("hits", "PropagatorAlg hits ntuple");
  if (m_ths->regTree("/rec/NtuplesHits", m_outputTree).isFailure()) {
    log << MSG::ERROR << "Couldn't register hits tree." << endmsg;
  }
  initializeTrees();

  log << MSG::INFO << "initialize fatrasSimAlg." << endmsg;

  return StatusCode::SUCCESS;
}
StatusCode FatrasSimulation::execute() {
  MsgStream log(msgSvc(), name());
  cleanTrees();

  auto tGeometry = m_geoSvc->trackingGeometry();

  using Navigator         = Acts::Navigator;
  using MagneticField     = Acts::ConstantBField;
  using ChargedStepper    = Acts::EigenStepper<>;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Navigator>;
  using NeutralStepper    = Acts::StraightLineStepper;
  using NeutralPropagator = Acts::Propagator<NeutralStepper, Navigator>;

  using Generator = std::mt19937;

  using CutPMin         = ActsFatras::Min<ActsFatras::Casts::P>;
  using ChargedSelector = CutPMin;

  using ChargedSimulation =
      ActsFatras::SingleParticleSimulation<ChargedPropagator, ActsFatras::StandardChargedElectroMagneticInteractions,
                                           HitSurfaceSelector, ActsFatras::NoDecay>;

  using NeutralSelector     = CutPMin;
  using NeutralInteractions = ActsFatras::InteractionList<ActsFatras::PhotonConversion>;
  using NeutralSimulation   = ActsFatras::SingleParticleSimulation<NeutralPropagator, NeutralInteractions,
                                                                 ActsFatras::NoSurface, ActsFatras::NoDecay>;

  using Simulation = ActsFatras::Simulation<ChargedSelector, ChargedSimulation, NeutralSelector, NeutralSimulation>;

  Acts::GeometryContext      geoCtx;
  Acts::MagneticFieldContext magCtx;
  Acts::Logging::Level       logLevel = Acts::Logging::Level::DEBUG;

  auto trackingGeometry = m_geoSvc->trackingGeometry();

  Navigator      navigator({trackingGeometry});
  ChargedStepper chargedStepper(
      std::make_shared<Acts::ConstantBField>(Acts::Vector3{0, 0, 2 * Acts::UnitConstants::T}));
  ChargedPropagator chargedPropagator(std::move(chargedStepper), navigator);
  NeutralPropagator neutralPropagator(NeutralStepper(), navigator);

  ChargedSimulation simulatorCharged(std::move(chargedPropagator),
                                     Acts::getDefaultLogger("ChargedSimulation", logLevel));
  NeutralSimulation simulatorNeutral(std::move(neutralPropagator),
                                     Acts::getDefaultLogger("NeutralSimulation", logLevel));
  Simulation        simulator(std::move(simulatorCharged), std::move(simulatorNeutral));

  Generator generator;

  std::vector<ActsFatras::Particle> simulatedInitial;
  std::vector<ActsFatras::Particle> simulatedFinal;
  std::vector<ActsFatras::Hit>      hits;

  const SimParticleContainer* p_partvect     = p_partvec.get();
  const auto&                 inputParticles = p_partvect;

  auto result = simulator.simulate(geoCtx, magCtx, generator, *inputParticles, simulatedInitial, simulatedFinal, hits);

  for (auto& hit : hits) {
    m_x.push_back(hit.position().x());
    m_y.push_back(hit.position().y());
    m_z.push_back(hit.position().z());
  }

  m_outputTree->Fill();

  return StatusCode::SUCCESS;
}
StatusCode FatrasSimulation::finalize() { return StatusCode::SUCCESS; }

StatusCode FatrasSimulation::initializeTrees() {
  m_outputTree->Branch("g_x", &m_x);
  m_outputTree->Branch("g_y", &m_y);
  m_outputTree->Branch("g_z", &m_z);

  return StatusCode::SUCCESS;
}

StatusCode FatrasSimulation::cleanTrees() {
  m_x.clear();
  m_y.clear();
  m_z.clear();

  return StatusCode::SUCCESS;
}
