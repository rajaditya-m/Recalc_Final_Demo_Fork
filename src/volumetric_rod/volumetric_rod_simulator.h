#ifndef VOLUMETRIC_ROD_SIMULATOR_H
#define VOLUMETRIC_ROD_SIMULATOR_H

#include "input_handler.h"
#include "singleton.h"

class VolumetricRodSimulator : public InputHandler
{
public:
  virtual int Idle();
  virtual int Render();
  virtual int HandleKeyPress(QKeyEvent *e);

private:
  VolumetricRodSimulator();
  class SimulatorInternal;
  SimulatorInternal* simulator_;
  virtual ~VolumetricRodSimulator();
  DECLARE_SINGLETON_CLASS(VolumetricRodSimulator)
};

#endif // VOLUMETRIC_ROD_SIMULATOR_H
