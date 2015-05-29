#include <memory>
#include "pbd_rod_simulator.h"
#include "pbd_elastic_rod.h"

class PBDRodSimulator::SimulatorInternal
{
public:
  SimulatorInternal()
    : rod_(new PBDElasticRod(rod::rod_v_num, rod::rod_length))
  {
  }

  ~SimulatorInternal()
  {
    delete rod_;
  }
  PBDElasticRod* rod_;
};

PBDRodSimulator::PBDRodSimulator()
  : InputHandler(0)
  , simulator_(new SimulatorInternal)
{
}

int PBDRodSimulator::Idle()
{
  for (int i = 0; i < global::simulation_step_per_idle; ++i) {
    simulator_->rod_->Simulate(global::time_step);
  }
  return -1;
}

int PBDRodSimulator::Render()
{
  simulator_->rod_->Render();
  return -1;
}

PBDRodSimulator::~PBDRodSimulator()
{
  delete simulator_;
}
