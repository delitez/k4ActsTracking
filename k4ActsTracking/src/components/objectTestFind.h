#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"

class objectTestFind : public GaudiAlgorithm {
public:
  explicit objectTestFind(const std::string&, ISvcLocator*);
  virtual ~objectTestFind();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:

};

struct MyTestVectorS  : public DataObject {
  std::vector<double> m_MyTestVector;
};
