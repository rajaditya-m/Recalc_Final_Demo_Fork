#ifndef SIMULATOR_H
#define SIMULATOR_H
#include "input_handler.h"

class PBDRodSimulator : public InputHandler
{
public:
  PBDRodSimulator();
  virtual ~PBDRodSimulator();
  virtual int Idle();
  virtual int Render();
  class SimulatorInternal;
  SimulatorInternal* simulator_;
  DECLARE_SINGLETON_CLASS(PBDRodSimulator)
};

#endif // SIMULATOR_H
