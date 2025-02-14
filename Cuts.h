#ifndef CUTS_H
#define CUTS_H

#include "TEnv.h"

class Cuts()
{
 public:

  Cuts(const std::string configfile)
    {
      m_penv = new TEnv(configfile.c_str());
    };
  
  ~Cuts();

  void setConfigFile(const std::string configfile)
  {
    m_penv->ReadFile(configfile);
  }


 private:

  // Configuration
  TEnv *m_penv{nullptr};

  
  
  
}

#endif
