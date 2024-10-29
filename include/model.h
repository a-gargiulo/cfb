#ifndef MODEL_H
#define MODEL_H

#include "parsing.h"
#include <string>
#include <vector>
#include <unordered_map>

namespace model
{

extern const double atm2pa;

struct MoleFraction
{
    std::string element;
    double fraction;

    MoleFraction(std::string elm, double frac) : element(elm), fraction(frac)
    {
    }
};

struct ThermoState
{
    double p, T;
    std::vector<MoleFraction> X;
};

class CounterflowBurner
{
  public:
    virtual void compute_exhaust_thermo_state(const parsing::InputData& inputs, ThermoState& state, int& err) = 0;
    
  protected:
    static void process_composition_string(const std::string& str, std::vector<MoleFraction>& compStr); 
};

class Equilibrium : public CounterflowBurner
{
  public:
    Equilibrium(const parsing::InputData& inputs);

    void compute_exhaust_thermo_state(const parsing::InputData& inputs, ThermoState& state, int& err) override;

  private:
    ThermoState inflow;
    ThermoState exhaust;

  private:
    void determine_reactants_composition(const parsing::InputData& inputs, std::unordered_map<std::string, double>& reacts, std::string& reactsStr, int& err) const;

};

} // namespace model

#endif // MODEL_H
