#include "logging.h"
#include "model.h"
#include "parsing.h"
#include <cantera/core.h>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>
#include <vector>


namespace model
{

Equilibrium::Equilibrium(const parsing::InputData& inputs) : CounterflowBurner(inputs)
{
}

void Equilibrium::compute_exhaust_thermo_state(const parsing::InputData& inputs, int& err)
{
    std::unordered_map<std::string, double> reacts;
    reacts.reserve(2 * NUM_EXPECTED_CORE_ELEMENTS);
    std::string reactsStr;

    std::unordered_map<std::string, double> prods;
    prods.reserve(100);

    determine_reactants_composition(inputs, reacts, reactsStr, err);
    return;
}

void Equilibrium::determine_reactants_composition(const parsing::InputData& inputs,
                                                  std::unordered_map<std::string, double>& reacts,
                                                  std::string& reactsStr, int& err) const
{
    Logger& logger = Logger::get_instance();

    std::vector<MoleFraction> fCore;
    fCore.reserve(NUM_EXPECTED_CORE_ELEMENTS);
    std::unordered_map<std::string, int> fCHO;
    fCHO.reserve(3);

    std::vector<MoleFraction> oxCore;
    oxCore.reserve(NUM_EXPECTED_CORE_ELEMENTS);

    std::string fCoreStr = std::get<std::string>(inputs.at("fuel_core_flow_composition"));
    std::string oxCoreStr = std::get<std::string>(inputs.at("oxidizer_core_flow_composition"));

    process_composition_string(fCoreStr, fCore);
    process_composition_string(oxCoreStr, oxCore);

    if (oxCore.begin()->element != "O2")    
    {
        std::ostringstream errMsg;
        errMsg << "The first species in the provided input oxidizer core flow composition is "
               << oxCore.begin()->element
               << ". The string must start with the oxidizer, O2, followd by all diluent species.";
        logger.log(LogLevel::ERROR, errMsg.str());
        err = -1;
        return;
    }

    



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
