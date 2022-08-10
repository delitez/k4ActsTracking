#include "EventGeneratorAlg.h"
#include "GaudiKernel/Service.h"


using namespace Gaudi;

DECLARE_COMPONENT(EventGeneratorAlg)

EventGeneratorAlg::EventGeneratorAlg(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

EventGeneratorAlg::~EventGeneratorAlg() { }

StatusCode EventGeneratorAlg::initialize() {

  return StatusCode::SUCCESS;
}
StatusCode EventGeneratorAlg::execute() {

   SimParticleContainer particles;

   static int seed = 5;
   std::mt19937 rng{seed++};
   std::normal_distribution<double> gauss(0., 1.);
   std::round(gauss(rng));

  size_t nPrimaryVertices = 0;

  for (size_t n = nMultiplicity; 0 < n; --n) {
    nPrimaryVertices+=1;

// TODO : Random 4vec does not work. It is set to 0,0,0,0 temporarily.

//    auto vertexPosition = (*vertex)(rng);
Acts::Vector4 vertexPosition(0., 0., 0., 0.);


    auto updateParticleInPlace = [&](ActsFatras::Particle& particle) {

        const auto pid = ActsFatras::Barcode(particle.particleId())
                                 .setVertexPrimary(nPrimaryVertices);

        const auto pos4 = (vertexPosition + particle.fourPosition()).eval();

        particle = particle.withParticleId(pid).setPosition4(pos4);
    };


     auto vertexParticles = genVertexParticles(rng,gauss);

    for (auto& vertexParticle : vertexParticles) {
      updateParticleInPlace(vertexParticle);
    }

    particles.merge(std::move(vertexParticles));

  }


  Container *containedParticle = new Container();
  containedParticle->m_particles = particles;


// This registiration must be revised because path is given in the python code and it crashes when it is not changed each time ---- TODO: Find a better way
// Question: Should particles registered individually or all together  --- if all together, then the EventGeneratorAlg must run only once and it is okay as it is now.
// Remark: Same name is okay if no compile, if compile, same name causes a crash

  StatusCode sc = eventSvc()->registerObject(objectPath, containedParticle);

  if( sc.isFailure() ) {
     std::cout << "Object cannot be registered." << std::endl;
     return StatusCode::FAILURE;
  }
  else{
    std::cout << "Object is registered in " << objectPath << std::endl;
  }

  return StatusCode::SUCCESS;
}
StatusCode EventGeneratorAlg::finalize() {

  return StatusCode::SUCCESS;
}


SimParticleContainer EventGeneratorAlg::genVertexParticles(std::mt19937& rng, std::normal_distribution<double>& gauss) {

  using UniformIndex = std::uniform_int_distribution<unsigned int>;


  std::uniform_real_distribution<double> phiDist(0, 2 * M_PI);
  std::uniform_real_distribution<double> etaDist(-4.0, 4.0);
  std::uniform_real_distribution<double> ptDist(10., 20.);
  std::uniform_real_distribution<double> qDist(0., 1.);

  double d0     = d0Sigma * gauss(rng);
  double z0     = z0Sigma * gauss(rng);
  double charge = qDist(rng) > 0.5 ? 1. : -1.;
  double t      = tSigma * gauss(rng);

  UniformIndex particleTypeChoice(0u, qDist(rng) ? 1u : 0u);
  Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;

  const Acts::PdgParticle pdgChoices[] = {pdg, static_cast<Acts::PdgParticle>(-pdg), };    //warum leer zeichen?
                                                                                          // braucht man 3 elements statt 2?
  const double qChoices[] = {charge,-charge, };

  SimParticleContainer particles;

  for (size_t ip = 1 ; ip <= nParticles; ++ip) {
    const auto pid = ActsFatras::Barcode(0u).setParticle(ip);
    const unsigned int type = particleTypeChoice(rng);
    const double q = qChoices[type];
    const double phi    = phiDist(rng);
    double eta    = etaDist(rng);
    double theta  = 2 * atan(exp(-eta));
    double pt     = ptDist(rng);
    double p      = pt / sin(theta);

    Acts::Vector3 direction;

    direction[Acts::eMom0] = sin(theta) * cos(phi);
    direction[Acts::eMom1] = sin(theta) * sin(phi);
    direction[Acts::eMom2] = cos(theta);

    //TODO: Fix the mass -- how to extract it from pdg? config file? by hand?

    double mass = 0.5;

    ActsFatras::Particle particle(pid, pdg, q, mass);
    particle.setDirection(direction);
    particle.setAbsoluteMomentum(p);

    particles.insert(particles.end(), std::move(particle));

  }

return particles;

};
