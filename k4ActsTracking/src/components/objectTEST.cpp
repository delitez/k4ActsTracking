#include "objectTEST.h"

#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>
#include "PropagatorAlg.h"

//using std::filesystem::current_path;


DECLARE_COMPONENT(objectTest)

objectTest::objectTest(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {}

objectTest::~objectTest() {}

StatusCode objectTest::initialize() {return StatusCode::SUCCESS; }

StatusCode objectTest::execute() {

  typedef std::vector<int> MyTestVector;
  DataObject *pObject;
  //std::string objectPath = "./testDir";

  StatusCode sc = eventSvc()->registerObject(".", pObject);

   if( sc.isFailure() ) {
      std::cout << "CANNOT initialize sc" << std::endl;
      return StatusCode::FAILURE;
   }
   else{
     std::cout << "CAN initialize sc" << std::endl;
   }

   MyTestVector *tv = 0;
   tv = dynamic_cast<MyTestVector *> (pObject);


   std::cout << "Object Test is alive!" << std::endl;
  return StatusCode::SUCCESS;
}

StatusCode objectTest::finalize() { return StatusCode::SUCCESS; }
