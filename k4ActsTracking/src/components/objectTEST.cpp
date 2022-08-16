#include "objectTEST.h"

#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>
#include "PropagatorAlg.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
//using std::filesystem::current_path;


DECLARE_COMPONENT(objectTest)

objectTest::objectTest(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

objectTest::~objectTest() {}

StatusCode objectTest::initialize() {return StatusCode::SUCCESS; }

StatusCode objectTest::execute() {

  typedef std::vector<int> MyTestVector;
  DataObject *pObject;
  MyTestVectorS *m_vector = new MyTestVectorS();



  StatusCode sc = eventSvc()->registerObject("/Event/Test", m_vector);
  if( sc.isFailure() ) {
     std::cout << "CANNOT initialize sc" << std::endl;
     return StatusCode::FAILURE;
  }
  else{
    std::cout << "CAN initialize sc" << std::endl;
  }


   std::cout << "Object Test is alive!" << std::endl;
  return StatusCode::SUCCESS;
}

StatusCode objectTest::finalize() { return StatusCode::SUCCESS; }
