// D. Elitez, August 2022
// Particle gun for gaudi4acts

#include "ParticleGunAlg.h"

using namespace Gaudi;

DECLARE_COMPONENT(ParticleGunAlg)

ParticleGunAlg::ParticleGunAlg(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

ParticleGunAlg::~ParticleGunAlg() {}

StatusCode ParticleGunAlg::initialize() {
  auto gen_gauss = randSvc()->generator(Rndm::Gauss(0., 1.));
  auto gen_flat  = randSvc()->generator(Rndm::Flat(0., 1.));

  return StatusCode::SUCCESS;
}
StatusCode ParticleGunAlg::execute() {
  MsgStream log(msgSvc(), name());

  size_t nPrimaryVertices = 0;

  for (size_t n = nMultiplicity; 0 < n; --n) {
    nPrimaryVertices += 1;

    // TODO : Random 4vec does not work. It is set to 0,0,0,0 temporarily.
    // Solution: set stuff per hand as rand gauss

    //    auto vertexPosition = (*vertex)(rng);
    auto vertexPosition  = Acts::Vector4(0., 0., 0., 0.);
    auto vertexParticles = genVertexParticles();

    Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;

    auto updateParticleInPlace = [&](ActsFatras::Particle& particle) {  // aus truth level zu vertexpos
      // std::cout << "particle before update: " << particle << std::endl;
      // std::cout << "particle 4 position before update: " << particle.fourPosition() << std::endl;
      // std::cout << "nPrimaryVertices: " << nPrimaryVertices << std::endl;
      //   const auto pid = ActsFatras::Barcode(particle.particleId()).setVertexPrimary(nPrimaryVertices);
      const auto pid = ActsFatras::Barcode(0u).setParticle(nPrimaryVertices);
      // std::cout << "Define pid " << std::endl;
      // std::cout << "vertexPosition: " << vertexPosition << std::endl;
      // std::cout << "particle four position: " << particle.fourPosition() << std::endl;
      // std::cout << "pos4 evaluation: " << (vertexPosition + particle.fourPosition()).eval() << std::endl;
      //
      // const auto pos4 = (vertexPosition + particle.fourPosition()).eval();
      // std::cout << "Define pos4 " << std::endl;
      // particle = particle.withParticleId(pid).setPosition4(pos4);
      particle.setParticleId(pid);
      particle.setPosition4(1., 1., 1., 1.);
      // std::cout << "particle after update: " << particle << std::endl;
      // std::cout << "particle 4 position after update: " << particle.fourPosition() << std::endl;
      // std::cout << "nPrimaryVertices: " << nPrimaryVertices << std::endl;

    };

    for (auto& vertexParticle : vertexParticles) {
      updateParticleInPlace(vertexParticle);
    }
    particles.merge(std::move(vertexParticles));
  }

  m_partvec.put(std::move(particles));

  return StatusCode::SUCCESS;
}

StatusCode ParticleGunAlg::finalize() { return StatusCode::SUCCESS; }

SimParticleContainer ParticleGunAlg::genVertexParticles() {
  MsgStream log(msgSvc(), name());

  SimParticleContainer gen_particles;

  Rndm::Numbers gauss(randSvc(), Rndm::Gauss(0., 1.));
  Rndm::Numbers phiDist(randSvc(), Rndm::Flat(0., 2 * M_PI));
  Rndm::Numbers etaDist(randSvc(), Rndm::Flat(-4., 4.));
  Rndm::Numbers ptDist(randSvc(), Rndm::Flat(10., 20.));
  Rndm::Numbers qDist(randSvc(), Rndm::Flat(0., 1.));
  Rndm::Numbers partCharge(randSvc(), Rndm::Flat(-1., 1.));

  double d0     = d0Sigma * gauss();
  double z0     = z0Sigma * gauss();
  double charge = qDist() > 0.5 ? 1. : -1.;
  double t      = tSigma * gauss();

  Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;

  const Acts::PdgParticle pdgChoices[] = {
      pdg,
      static_cast<Acts::PdgParticle>(-pdg),
  };
  const double qChoices[] = {
      charge,
      -charge,
  };

  for (size_t ip = 1; ip <= nParticles; ++ip) {
    const auto pid = ActsFatras::Barcode(0u).setParticle(ip);

    double q = 0.;

    if (charge > 0) {
      q   = 1.;
      pdg = pdg;
    } else if (charge < 0) {
      q   = -1.;
      pdg = static_cast<Acts::PdgParticle>(-pdg);
    } else {
      log << MSG::ERROR << "Cannot set q and pdg." << endmsg;
      StatusCode::FAILURE;
    }

    const double phi   = phiDist();
    double       eta   = etaDist();
    double       theta = 2 * atan(exp(-eta));
    double       pt    = ptDist();
    double       p     = pt / sin(theta);

    Acts::Vector3 direction;

    direction[Acts::eMom0] = sin(theta) * cos(phi);
    direction[Acts::eMom1] = sin(theta) * sin(phi);
    direction[Acts::eMom2] = cos(theta);

    double mass = ActsFatras::findMass(pdg);

    ActsFatras::Particle _particle(pid, pdg, q, mass);
    _particle.setDirection(direction);
    _particle.setAbsoluteMomentum(p);

    gen_particles.insert(gen_particles.end(), std::move(_particle));
  }

  return gen_particles;
};
