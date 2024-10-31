#ifndef MODEL_H
#define MODEL_H

#include "parsing.h"
#include <string>
#include <vector>
#include <unordered_map>

#define NUM_EXPECTED_CORE_ELEMENTS 5

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
    CounterflowBurner(const parsing::InputData& inputs);

    virtual void compute_exhaust_thermo_state(const parsing::InputData& inputs, int& err) = 0;
    
  protected:
    static void process_composition_string(const std::string& str, std::vector<MoleFraction>& compStr); 
    static bool fuel_is_hydrocarbon(const std::string& fuel, std::unordered_map<std::string, int>& nElements);

  protected:
    ThermoState inflow;
    ThermoState exhaust;
};

class Equilibrium : public CounterflowBurner
{
  public:
    Equilibrium(const parsing::InputData& inputs);

    void compute_exhaust_thermo_state(const parsing::InputData& inputs, int& err) override;

  private:
    void determine_reactants_composition(const parsing::InputData& inputs, std::unordered_map<std::string, double>& reacts, std::string& reactsStr, int& err) const;

};

} // namespace model

#endif // MODEL_H
