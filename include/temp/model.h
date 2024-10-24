#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>

#include "utils.h"

namespace model
{

struct MoleFraction
{
    std::string element;
    double fraction;
}

struct ThermoState
{
    double p, T;
    std::vector<MoleFraction> X;
}

class CounterflowBurner
{
  public:
    virtual ThermoState& compute_exhaust_thermo_state(const std::vector<utils::KeyValue>& inputs) = 0;
}

class Equilibrium : public CounterflowBurner
{
  public:
    virtual ThermoState& compute_exhaust_thermo_state(const std::vector<utils::KeyValue>& inputs) override;

}

} // namespace model

#endif // MODEL_H
