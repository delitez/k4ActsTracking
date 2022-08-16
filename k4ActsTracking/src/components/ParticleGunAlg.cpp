#include "ParticleGunAlg.h"
#include "GaudiKernel/Service.h"


using namespace Gaudi;

DECLARE_COMPONENT(ParticleGunAlg)

ParticleGunAlg::ParticleGunAlg(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

ParticleGunAlg::~ParticleGunAlg() { }

StatusCode ParticleGunAlg::initialize() {

  return StatusCode::SUCCESS;
}
StatusCode ParticleGunAlg::execute() {



   static int seed = 5;
   std::mt19937 rng{seed++};
   std::normal_distribution<double> gauss(0., 1.);
   std::round(gauss(rng));

  size_t nPrimaryVertices = 0;

  for (size_t n = nMultiplicity; 0 < n; --n) {
    nPrimaryVertices+=1;

// TODO : Random 4vec does not work. It is set to 0,0,0,0 temporarily.
// Solution: set stuff per hand as rand gauss

//    auto vertexPosition = (*vertex)(rng);
     auto vertexPosition = Acts::Vector4(1.,1.,1.,1.);
     auto vertexParticles = genVertexParticles(rng,gauss);


Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;

     auto updateParticleInPlace = [&](ActsFatras::Particle& particle) { // aus truth level zu vertexpos
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
             particle.setPosition4(1.,1.,1.,1.);
             // std::cout << "particle after update: " << particle << std::endl;
             // std::cout << "particle 4 position after update: " << particle.fourPosition() << std::endl;
             // std::cout << "nPrimaryVertices: " << nPrimaryVertices << std::endl;

           };



           for (auto& vertexParticle : vertexParticles) {
             updateParticleInPlace(vertexParticle);
           }
   particles.merge(std::move(vertexParticles));
  }



    // std::cout << "BEFORE PUT WELT!" << std::endl;
     m_partvec.put(std::move(particles));
    // std::cout << "AFTER PUT WELT!" << std::endl;
    //
    //
    //




  return StatusCode::SUCCESS;
}


StatusCode ParticleGunAlg::finalize() {

  return StatusCode::SUCCESS;
}


SimParticleContainer ParticleGunAlg::genVertexParticles(std::mt19937& rng, std::normal_distribution<double>& gauss) {

  SimParticleContainer gen_particles;

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

    ActsFatras::Particle _particle(pid, pdg, q, mass);
    _particle.setDirection(direction);
    _particle.setAbsoluteMomentum(p);

    gen_particles.insert(gen_particles.end(), std::move(_particle));

  }

return gen_particles;

};
