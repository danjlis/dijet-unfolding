
#ifndef DIJETUNFOLD_H
#define DIJETUNFOLD_H

#include "TEnv.h"

class Cuts()
{
 public:
  DijetUnfold();
  ~DijetUnfold();

  int MakeResponseMatrix();

  int ClosureTest();

  int UnfoldData();

  void setConfigFile(const std::string configfile)
  {
    m_config_file = configfile;
  }

 private:

  // Configuration
  std::string m_config_file{""};

  // class to store my cuts
  Cuts *m_cuts{nullptr};
  
  
}

#endif
