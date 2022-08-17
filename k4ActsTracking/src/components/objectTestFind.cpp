#include "objectTestFind.h"

#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>
#include "PropagatorAlg.h"


//using std::filesystem::current_path;


DECLARE_COMPONENT(objectTestFind)

objectTestFind::objectTestFind(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

objectTestFind::~objectTestFind() {}

StatusCode objectTestFind::initialize() {return StatusCode::SUCCESS; }

StatusCode objectTestFind::execute() {

  
  typedef std::vector<int> MyTestVector;
  DataObject *pObject;
  MyTestVectorS *m_vector = new MyTestVectorS();


  StatusCode sc;

  sc = eventSvc()->retrieveObject("/Event/Test", pObject);
  if( sc.isFailure() ) {
     std::cout << "(no found) initialize sc" << std::endl;
    // return StatusCode::FAILURE;
    return sc;
  }
  else{
    std::cout << " (found) initialize sc" << std::endl;
  }



   std::cout << "Object Test Find is alive!" << std::endl;
  return StatusCode::SUCCESS;
}

StatusCode objectTestFind::finalize() { return StatusCode::SUCCESS; }
