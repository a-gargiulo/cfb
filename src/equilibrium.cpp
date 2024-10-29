#include "logging.h"
#include "model.h"
#include "parsing.h"
#include <cantera/core.h>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

namespace model
{

Equilibrium::Equilibrium(const parsing::InputData& inputs)
{
    auto sol = Cantera::newSolution(std::get<std::string>(inputs.at("reaction_mechanism")));
    auto gas = sol->thermo();

    size_t num_species = gas->nSpecies();

    inflow.p = std::get<double>(inputs.at("reactants_pressure"));
    inflow.T = std::get<double>(inputs.at("reactants_temperature"));
    inflow.X.reserve(num_species);

    for (size_t i = 0; i < num_species; ++i)
    {
        inflow.X.emplace_back(gas->speciesName(i), 0.0);
    }

    exhaust.p = std::get<double>(inputs.at("pressure"));
    exhaust.T = 0.0;
    exhaust.X.reserve(num_species);

    for (size_t i = 0; i < num_species; ++i)
    {
        exhaust.X.emplace_back(gas->speciesName(i), 0.0);
    }
}

void Equilibrium::compute_exhaust_thermo_state(const parsing::InputData& inputs, ThermoState& state, int& err)
{
    return;
}

void Equilibrium::determine_reactants_composition(const parsing::InputData& inputs,
                                                  std::unordered_map<std::string, double>& reacts,
                                                  std::string& reactsStr, int& err) const
{
    Logger& logger = Logger::get_instance();

    std::vector<MoleFraction> fCore;
    std::vector<MoleFraction> oxCore;

    std::string fCoreStr = std::get<std::string>(inputs.at("fuel_core_flow_composition"));
    std::string oxCoreStr = std::get<std::string>(inputs.at("oxidizer_core_flow_composition"));



    // if (oxCoreStr.begin()->first != "O2")
    // {
    //     std::ostringstream errMsg;
    //     errMsg << "The first species in the provided input oxidizer core flow composition is "
    //            << oxCoreStr.begin()->first
    //            << ". The string must start with the oxidizer, O2, followd by all diluent species.";
    //     logger.log(LogLevel::ERROR, errMsg.str());
    // }

    return;
}

} // namespace model
