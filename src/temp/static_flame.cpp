#define _USE_MATH_DEFINES #include<math.h>
#include "model.h"
#include "cantera/core.h"
#include "error.h"
#include "logging.h"
#include "strutil.h"
#include "utils.h"
#include <algorithm>
#include <any>
#include <cctype>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#define CHO_MAP_SIZE 3

namespace model
{

ErrorCode Cfb::compute_exhaust_temperature(const std::vector<utils::KeyAnyValue>& inputs)
{
    ErrorCode err;

    // Output variables
    std::unordered_map<std::string, double> reacts, prods;
    std::string reactsStr;
    std::unordered_map<std::string, double> reacts_adj, prods_adj;
    std::vector<double> ndot, mdot;

    // Local variables
    std::ostringstream oss;

    err = determine_reactants_composition(inputs, reacts, reactsStr);
    if (err != ErrorCode::SUCCESS) return err;

    err = compute_eq_products_composition(inputs, reactsStr, prods);
    if (err != ErrorCode::SUCCESS) return err;

    // err = flow_rates_from_mass_balance(inputs, eqProds, ndot, mdot);
    // if (err != ErrorCode::SUCCESS) return err;

    // err = adjust_reactants_mole_fractions(inputs, ndot, eqReactsComp, reacts_adj);
    // if (err != ErrorCode::SUCCESS) return err;

    // err = adjust_products_mole_fractions(inputs, ndot, eqProds, prods_adj);
    // if (err != ErrorCode::SUCCESS) return err;

    // // iterate
    // double qloss;
    // double Texhaust_old, Texhaust;
    // double Tguess = utils::get_value<double>(inputs, "initial_guess_exhaust_temperature");
    // double reltol = utils::get_value<double>(inputs, "relative_tolerance");
    // int maxIter = static_cast<int>(utils::get_value<double>(inputs, "maximum_number_of_iterations"));

    // Texhaust = Tguess;
    // Texhaust_old = 0.0;
    // int iter = 0;
    // while (std::fabs(Texhaust - Texhaust_old)/std::fabs(Texhaust_old) > reltol && iter < maxIter)
    // {
    //     Texhaust_old = Texhaust;

    //     // qloss = get_qloss_from_energy_balance(inputs, reacts_adj, prods_adj, nTot, Texhaust_old);
    //     qloss = get_qloss_from_heat_transfer(inputs, reacts_adj, prods_adj, nTot, Texhaust_old);

    //     // Texhaust = get_Texhaust_from_heat_transfer(inputs, reacts_adj, prods_adj, nTot, Texhaust_old, qloss);
    //     Texhaust = get_Texhaust_from_energy_balance(inputs, reacts_adj, prods_adj, nTot, mtot, qloss);
    //     iter++;
    //     std::cout<< Texhaust << std::endl;
    // }

    // double rho_exhaust = get_exhaust_density(inputs, prods_adj, Texhaust, nTot);
    // std::cout << "DENSITY: " << rho_exhaust << std::endl;
    // return Texhaust;
}

// double Cfb::get_exhaust_density(const std::vector<utils::KeyAnyValue>& inputs, const std::vector<utils::KeyAnyValue>&
// prods, double& Texhaust, const std::vector<double>& nTot)
// {
//    double rho_exhaust;

//    std::string rxn = utils::get_value<std::string>(inputs, "reaction_mechanism");
//    double p = utils::get_value<double>(inputs, "pressure");

//     auto sol = Cantera::newSolution(rxn);
//     auto gas = sol->thermo();

//     std::ostringstream pComp;
//     double nTotProds = nTot[2] + nTot[3] + nTot[4];
//     for (size_t i = 0; i < prods.size(); ++i)
//     {
//         pComp << prods[i].key << ": "  << std::any_cast<double>(prods[i].value) * nTotProds;
//         if (i == prods.size() - 1)
//             continue;
//         else
//             pComp << ", ";
//     }
//     gas->setState_TPX(Texhaust, p * Cantera::OneAtm, pComp.str());

//     rho_exhaust = gas->density();

//     return rho_exhaust;

// }

// double Cfb::get_qloss_from_heat_transfer(const std::vector<utils::KeyAnyValue>& inputs, const
// std::vector<utils::KeyAnyValue>& reacts, const std::vector<utils::KeyAnyValue>& prods, const std::vector<double>&
// nTot, double Tguess)
// {

//     std::string rxn = utils::get_value<std::string>(inputs, "reaction_mechanism");
//     double p = utils::get_value<double>(inputs, "pressure");
//     double Tinf = utils::get_value<double>(inputs, "ambient_temperature");

//     double Lc = utils::get_value<double>(inputs, "wall_length");
//     double t = utils::get_value<double>(inputs, "wall_thickness");

//     auto sol = Cantera::newSolution(rxn);
//     auto gas = sol->thermo();
//     auto trans = sol->transport();

//     // Nusselt Chamber
//     std::ostringstream pComp;
//     double nTotProds = nTot[2] + nTot[3] + nTot[4];
//     for (size_t i = 0; i < prods.size(); ++i)
//     {
//         pComp << prods[i].key << ": "  << std::any_cast<double>(prods[i].value) * nTotProds;
//         if (i == prods.size() - 1)
//             continue;
//         else
//             pComp << ", ";
//     }
//     gas->setState_TPX(Tguess, p * Cantera::OneAtm, pComp.str());
//     double kprods = trans->thermalConductivity();
//     double muprods = trans->viscosity();
//     double cpprods = gas->cp_mass();

//     double rhoprods = gas->density();
//     double Prprods = cpprods*muprods/kprods;

//     double vprods = utils::get_value<double>(inputs, "exhaust_gas_convective_velocity_magnitude");
//     double Reprods = rhoprods * vprods * Lc / muprods;

//     double h_prods = (kprods/Lc) * 0.664 * std::sqrt(Reprods) * std::pow(Prprods, 1.0/3.0);

//     // Ambient  Air
//     gas->setState_TPX(Tinf, 1 * Cantera::OneAtm, "O2:0.21, N2:0.78, AR:0.01");
//     double vAir = utils::get_value<double>(inputs, "ambient_air_convective_velocity_magnitude");
//     double kair = trans->thermalConductivity();
//     double muair = trans->viscosity();
//     double cpair = gas->cp_mass();
//     double rhoair = gas->density();
//     double Prair = cpair*muair/kair;
//     double Reair = rhoair * vAir * Lc / muair;

//     double h_air = (kair/Lc) * 0.664 * std::sqrt(Reair) * std::pow(Prair, 1.0/3.0);

//     // Wall Temperatures
//     double kwall = 45;
//     double Tinfw = h_air * (1-t*h_prods/kwall) * Tinf + h_prods * Tguess /(h_prods * (1-t*h_air/kwall) + h_air);

//     double Tcw = Tinfw*(1-t*h_air/kwall) + t*h_air*Tinf/kwall;

//     double qloss = -(Tcw - Tinfw) * kwall * Lc * Lc / t;
//     return qloss;

// }

// ErrorCode Cfb::adjust_reactants_mole_fractions(const std::vector<utils::KeyAnyValue>& inputs, const
// std::vector<double>& ndot, const std::string& rComp, std::vector<utils::KeyAnyValue>& reacts_adj)
// {
//     ErrorCode err,

//     std::vector<std::pair<std::string, double>> rCompVec = strutil::process_composition_string(rComp);
//     std::vector<utils::KeyAnyValue> reacts;
//     std::string fShroud;
//     std::string oxShroud;

//     [fShroud, err] = utils::get_value<std::string>(inputs, "fuel_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     [oxShroud, err] = utils::get_value<std::string>(inputs, "oxidizer_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     // Compute molar flow rate for 1 mol/s fuel
//     double rNdotTotNorm = 0 ;
//     for (size_t i = 0; i < rCompVec.size(); ++i)
//     {
//         rNdotTotNorm += rCompVec[i].second;
//     }

//     // Compute reactants mole fractions
//     for (size_t i = 0; i < rCompVec.size(); ++i)
//     {
//         reacts.push_back({rCompVec[i].first, rCompVec[i].second / rNdotTotNorm});
//     }

//     // Compute 'actual' molar flow rates
//     std::vector<double> rNdot;
//     for (size_t i = 0; i < reacts.size(); ++i)
//     {
//         rNdot.push_back((ndot[0] + ndot[1]) * std::any_cast<double>(reacts[i].value));
//     }

//     // Include shroud
//     size_t fShroudIdx = -1;
//     for (size_t i = 0; i < reacts.size(); ++i)
//     {
//         if (reacts[i].key == fShroud)
//         {
//             fShroudIdx = i;
//             rNdot[i] += ndot[2];
//             break;
//         }
//     }
//     if (idxFShroud == -1)
//     {
//         reacts.push_back({fShroud, 1.0}); // 1.0 is a dummy value
//         rNdot.push_back(ndot[2]);
//     }

//     size_t oxShroudIdx = -1;
//     for (size_t i = 0; i < reacts.size(); ++i)
//     {
//         if (reacts[i].key == oxShroud)
//         {
//             oxShroudIdx = i;
//             rNdot[i] += ndot[3];
//             break;
//         }
//     }
//     if (oxShroudIdx == -1)
//     {
//         reacts.push_back({oxShroud, 1.0}); // 1.0 is a dummy value
//         rNdot.push_back(ndot[3]);
//     }

//     for (size_t i = 0; i < reacts.size(); ++i)
//         reacts_adj.push_back({reacts[i].key, rNdot[i] / (ndot[0] + ndot[1] + ndot[2] + ndot[3])});

//     return ErrorCode::SUCCESS;
// }

// ErrorCode Cfb::adjust_products_mole_fractions(const std::vector<utils::KeyAnyValue>& inputs, const
// std::vector<double>& ndot, const std::vector<utils::KeyAnyValue>& prods, std::vector<utils::KeyAnyValue>& prods_adj)
// {
//     logging::Logger& logger = logging::Logger::get_instance();
//     std::string errorMessage;

//     std::vector<double> pNdot;
//     std::string fShroudSearch, oxShroudSearch;

//     double prodsNoShroudNdot = ndot[4];

//     [fShroudSearch, err] = utils::get_value<std::string>(inputs, "fuel_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     [oxShroudSearch, err] = utils::get_value<std::string>(inputs, "oxidizer_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     for (const utils::KeyAnyValue& kv : prods)
//     {
//         pNdot.push_back(prodsNoShroudNdot * std::any_cast<double>(kv.value));
//     }

//     // Add shroud gas influence to products and adjust mass fractions
//     size_t fShroudIdx = -1;
//     for (size_t i = 0; i < prods.size(); ++i) {
//         if (prods[i].key == fShroudSearch) {
//             fShroudIdx = i;
//             break;
//         }
//     }
//     if (fShroudIdx == -1)
//     {
//         errorMessage = (std::ostringstream() << "Shroud gas " << fShroud << " not found in the reaction file. Choose
//         a valid shroud gas.").str(); logger.log(logging::LogLevel::ERROR, errorMessage); return
//         ErrorCode::INPUT_INVALID;
//     }
//     else
//     {
//         pNdot[fShroudIdx] += ndot[2];
//     }

//     size_t oxShroudIdx = -1;
//     for (size_t i = 0; i < prods.size(); ++i) {
//         if (prods[i].key == oxShroudSearch) {
//             oxShroudIdx = i;
//             break;
//         }
//     }
//     if (oxShroudIdx == -1)
//     {

//         errorMessage = (std::ostringstream() << "Shroud gas " << fShroud << " not found in the reaction file. Choose
//         a valid shroud gas.").str(); logger.log(logging::LogLevel::ERROR, errorMessage); return
//         ErrorCode::INPUT_INVALID;
//     }
//     else
//     {
//         pNdot[oxShroudIdx] += ndot[3];
//     }

//     for (size_t i = 0; i < pNdot.size(); ++i)
//     {
//         prods_adj.push_back({prods[i].key, pNdot[i] / ndot[5]});
//     }

//     return ErrorCode::SUCCESS;
// }

bool Cfb::fuel_is_hydrocarbon(const std::string& fuel, std::unordered_map<std::string, int>& elementCounts)
{
    std::stack<std::unordered_map<std::string, int>> parenStack;

    size_t length = fuel.size();
    size_t i = 0;

    while (i < length)
    {
        char ch = fuel[i];

        if (ch == 'C' || ch == 'H' || ch == 'O')
        {
            // Parse element and optional count
            size_t j = i + 1;
            while (j < length && strutil::is_digit(fuel[j]))
            {
                j++;
            }

            std::string element(1, ch);
            std::string countStr = fuel.substr(i + 1, j - i - 1);
            int count = countStr.empty() ? 1 : std::stoi(countStr);

            if (count <= 0)
            {
                return false; // Invalid count
            }

            elementCounts[element] += count;
            i = j;
        }
        else if (ch == '(')
        {
            // Push current element counts onto the stack and start a new scope
            parenStack.push(elementCounts);
            elementCounts.clear();
            i++;
        }
        else if (ch == ')')
        {
            // Pop the stack and handle the multiplier
            if (parenStack.empty())
            {
                return false; // Unmatched parenthesis
            }

            std::unordered_map<std::string, int> innerElements = elementCounts;
            elementCounts = parenStack.top();
            parenStack.pop();

            i++;
            // Handle multiplier
            size_t j = i;
            while (j < length && strutil::is_digit(fuel[j]))
            {
                j++;
            }

            std::string multiplierStr = fuel.substr(i, j - i);
            int multiplier = multiplierStr.empty() ? 1 : std::stoi(multiplierStr);

            if (multiplier <= 0)
            {
                return false; // Invalid multiplier
            }

            // Multiply inner elements
            for (const auto& pair : innerElements)
            {
                elementCounts[pair.first] += pair.second * multiplier;
            }

            i = j;
        }
        else if (strutil::is_digit(ch))
        {
            // Invalid state: multipliers should directly follow a closing parenthesis, not standalone
            return false;
        }
        else
        {
            // Invalid character
            return false;
        }
    }

    // Check if all opened parentheses are closed
    return parenStack.empty();
}

ErrorCode Cfb::determine_reactants_composition(const std::vector<utils::KeyAnyValue>& inputs,
                                               std::unordered_map<std::string, double>& reacts,
                                               std::string& reactsStr)
{
    logging::Logger& logger = logging::Logger::get_instance();

    std::string errMsg;

    std::unordered_map<std::string, int> fCHO;
    std::unordered_map<std::string, double> rStoichCoeffs;

    std::vector<std::pair<std::string, double>> fCoreComp, oxCoreComp;

    auto [fCoreCompStr, err1] = utils::get_value<std::string>(inputs, "fuel_core_flow_composition");
    if (err1 != ErrorCode::SUCCESS)
        return err1;

    auto [oxCoreCompStr, err2] = utils::get_value<std::string>(inputs, "oxidizer_core_flow_composition");
    if (err2 != ErrorCode::SUCCESS)
        return err2;

    if (oxCoreComp.begin()->first != "O2")
    {
        errMsg =
            (std::ostringstream() << "The first species in the provided input oxidizer core flow composition is "
                                  << oxCoreComp.begin()->first
                                  << ". The string must start with the oxidizer, O2, followed by all diluent species.")
                .str();
        logger.log(logging::LogLevel::ERROR, errMsg);
        return ErrorCode::INPUT_INVALID;
    }

    fCoreComp = strutil::process_composition_string(fCoreCompStr);
    oxCoreComp = strutil::process_composition_string(oxCoreCompStr);

    if (!fuel_is_hydrocarbon(fCoreComp.begin()->first, fCHO))
    {
        errMsg =
            (std::ostringstream() << fCoreComp.begin()->first
                                  << " is not a hydrocarbon. Select a valid hydrocarbon fuel and ensure that it is the "
                                     "first species provided in the input fuel core flow composition string.")
                .str();
        logger.log(logging::LogLevel::ERROR, errMsg);
        return ErrorCode::INPUT_INVALID;
    }

    // TO-DO: Add the use of equivalence ratio
    rStoichCoeffs[fCoreComp.begin()->first] = 1;
    rStoichCoeffs[oxCoreComp.begin()->first] = fCHO["C"] + fCHO["H"] / 4.0 - fCHO["O"] / 2.0;

    // Handle diluents
    if (oxCoreComp.size() > 1)
    {
        auto oxCoreIt = oxCoreComp.begin();
        std::advance(oxCoreIt, 1);

        for (; oxCoreIt != oxCoreComp.end(); ++oxCoreIt)
        {
            rStoichCoeffs[oxCoreIt->first] =
                rStoichCoeffs[oxCoreComp.begin()->first] * oxCoreIt->second / oxCoreComp.begin()->second;
        }
    }

    if (fCoreComp.size() > 1)
    {
        auto fCoreIt = fCoreComp.begin();
        std::advance(fCoreIt, 1);

        std::string searchKey;
        std::vector<std::pair<std::string, double>>::iterator searchIt;

        for (; fCoreIt != fCoreComp.end(); ++fCoreIt)
        {
            searchKey = fCoreIt->first;
            searchIt =
                std::find_if(oxCoreComp.begin(), oxCoreComp.end(),
                             [&searchKey](const std::pair<std::string, double>& p) { return p.first == searchKey; });
            // Fuel diluent species is also an oxidizer diluent species
            if (searchIt != oxCoreComp.end())
            {
                rStoichCoeffs[fCoreIt->first] +=
                    rStoichCoeffs[fCoreComp.begin()->first] * fCoreIt->second / fCoreComp.begin()->second;
            }
            // Fuel diluent species has not been encountered yet
            else
            {
                rStoichCoeffs[fCoreIt->first] =
                    rStoichCoeffs[fCoreComp.begin()->first] * fCoreIt->second / fCoreComp.begin()->second;
            }
        }
    }

    double rNdotTot = 0;
    for (const auto& p : rStoichCoeffs)
        rNdotTot += p.second;

    for (const auto& p : rStoichCoeffs)
        reacts[p.first] = p.second / rNdotTot;

    std::ostringstream oss;
    for (const auto& p : rStoichCoeffs)
    {
        oss << p.first << ": " << p.second << ", ";
    }
    reactsStr = oss.str();
    reactsStr.erase(reactsStr.end() - 2, reactsStr.end()); 

    return ErrorCode::SUCCESS;
}

ErrorCode Cfb::compute_eq_products_composition(const std::vector<utils::KeyAnyValue>& inputs, const std::unordered_map<std::string, double>& reacts, std::unordered_map<std::string, double>& prods)
{
    auto [rxn, err1] = utils::get_value<std::string>(inputs, "reaction_mechanism");
    if (err1 != ErrorCode::SUCCESS) return err1;
    auto [p, err2] = utils::get_value<double>(inputs, "pressure");
    if (err2 != ErrorCode::SUCCESS) return err2;
    auto [Tr, err3] = utils::get_value<double>(inputs, "reactants_temperature");
    if (err3 != ErrorCode::SUCCESS) return err3;

    auto sol = Cantera::newSolution(rxn);
    auto gas = sol->thermo();

    gas->setState_TPX(Tr, p * Cantera::OneAtm, eqReactsComp);
    gas->equilibrate("HP");

    for (size_t i = 0; i < gas->nSpecies(); ++i)
    {
        eqProds.push_back({gas->speciesName(i), gas->moleFraction(i)});
    }

    return ErrorCode::SUCCESS;
}

// ErrorCode Cfb::flow_rates_from_mass_balance(const std::vector<utils::KeyAnyValue>& inputs, const
// std::vector<utils::KeyAnyValue>& eqProds, std::vector<double>& ndot, std::vector<double>& mdot)
// {
//     ErrorCode err;

//     size_t idx;

//     double aG;
//     double p, Tr;
//     double nozzleDi, nozzleL;
//     std::string rxn;
//     std::string fCoreCompStr, oxCoreCompStr;
//     std::string fShourd, oxShroud;
//     double momFrac;
//     std::vector<std::pair<std::string, double>> fCoreComp, oxCoreComp;
//     std::vector<double> fCoreM, oxCoreM, eqProdsM;

//     double nozzleAi = Di * Di * M_PI / 4;
//     double nozzleAo = M_PI * (std::pow(0.01495 / 2, 2) - std::pow(0.00726 / 2, 2));
//     double Ru = 8.314;

//     [aG, err] = utils::get_value<double>(inputs, "global_strain_rate");
//     if (err != ErrorCode::SUCCESS) return err;

//     [Tr, err] = utils::get_value<double>(inputs, "reactants_temperature");
//     if (err != ErrorCode::SUCCESS) return err;

//     [p, err] = utils::get_value<double>(inputs, "pressure");
//     if (err != ErrorCode::SUCCESS) return err;

//     [nozzleDi, err] = utils::get_value<double>(inputs, "inner_nozzle_diameter");
//     if (err != ErrorCode::SUCCESS) return err;

//     [nozzleL, err] = utils::get_value<double>(inputs, "nozzle_separation_distance");
//     if (err != ErrorCode::SUCCESS) return err;

//     [fCoreCompStr, err] = utils::get_value<std::string>(inputs, "fuel_core_flow_composition");
//     if (err != ErrorCode::SUCCESS) return err;

//     [oxCoreCompStr, err] = utils::get_value<std::string>(inputs, "oxidizer_core_flow_composition");
//     if (err != ErrorCode::SUCCESS) return err;

//     [fShroud, err] = utils::get_value<std::string>(inputs, "fuel_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     [oxShroud, err] = utils::get_value<std::string>(inputs, "oxidizer_shroud");
//     if (err != ErrorCode::SUCCESS) return err;

//     [momFrac, err] = utils::get_value<double>(inputs, "shroud_momentum_fraction");
//     if (err != ErrorCode::SUCCESS) return err;

//     [rxn, err] = utils::get_value<std::string>(inputs, "reaction_mechanism");
//     if (err != ErrorCode::SUCCESS) return err;

//     fCoreComp = strutil::process_composition_string(fCompStr);
//     oxCoreComp = strutil::process_composition_string(oxCompStr);

//     auto sol = Cantera::newSolution(rxn);
//     auto gas = sol->thermo();

//     // Retrieve molecular weights
//     for (const std::pair<std::string, double>& elm : fCoreComp)
//     {
//        idx = gas->speciesIndex(elm.first);
//        fCoreM.push_back(gas->molecularWeight(idx) / 1000.0);
//     }

//     for (const std::pair<std::string, double>& elm : oxCoreComp)
//     {
//        idx = gas->speciesIndex(elm.first);
//        oxCoreM.push_back(gas->molecularWeight(idx) / 1000.0);
//     }

//     for (const utils::KeyAnyValue& elm : eqProds)
//     {
//        idx = gas->speciesIndex(elm.key);
//        eqProdsM.push_back(gas->molecularWeight(idx) / 1000.0);
//     }

//     idx = gas->speciesIndex(fShroud);
//     double fShroudM = gas->molecularWeight(idx) / 1000.0;

//     idx = gas->speciesIndex(oxShroud);
//     double oxShroudM = gas->molecularWeight(idx) / 1000.0;

//     double fCoreMavg = 0;
//     std::vector<double>::iterator fCoreMit = fCoreM.begin();
//     for (const std::pair<std::string, double>& elm : fCoreComp)
//     {
//         fCoreMavg += elm.second * *fCoreMit;
//         fCoreMit++;
//     }

//     double oxCoreMavg = 0;
//     std::vector<double>::iterator oxCoreMit = oxCoreM.begin();
//     for (const std::pair<std::string, double>& elm : oxCoreComp)
//     {
//         oxCoreMavg += elm.second * *oxCoreMit;
//         oxCoreMit++;
//     }

//     double eqProdsMavg = 0;
//     std::vector<double>::iterator eqProdsMit = eqProdsM.begin();
//     for (const utils::KeyAnyValue& kv : eqProds)
//     {
//         eqProdsMavg += std::any_cast<double>(kv.value) * *eqProdsMit;
//         eqProdsMit++;
//     }

//     // Compute molar and mass flow rates
//     double oxCoreV = nozzleL * aG / 4;
//     double oxCoreVdot = oxCoreV * nozzleAi;
//     double oxCoreNdot = p * Cantera::OneAtm * oxCoreVdot / Ru / Tr;
//     double oxCoreMdot = oxCoreNdot * oxCoreMavg;

//     double fCoreRho = p * Cantera::OneAtm * fCoreMavg / Ru / Tr;
//     double oxCoreRho = p * Cantera::OneAtm * oxCoreMavg / Ru / Tr;

//     double fCoreV = oxCoreV * std::sqrt(oxCoreRho / fCoreRho);
//     double fCoreVdot = fCoreV * nozzleAi;
//     double fCoreNdot = p * Cantera::OneAtm * fCoreVdot / Ru / Tr;
//     double fCoreMdot = fCoreNdot * fCoreMavg;

//     double fShroudRho = p * Cantera::OneAtm * fShroudM / Ru / Tr;
//     double oxShroudRho = p * Cantera::OneAtm * oxShroudM / Ru / Tr;

//     double oxShroudV = oxCoreV * std::sqrt(momFrac * oxCoreRho / oxShroudRho);
//     double oxShroudVdot = oxShroudV * nozzleAo;
//     double oxShroudNdot = p * Cantera::OneAtm * oxShroudVdot / Ru / Tr;
//     double oxShroudMdot = oxShroudNdot * oxShroudM;

//     double fShroudV = oxShroudV * std::sqrt(oxShroudRho / fShroudRho);
//     double fShroudVdot = fShroudV * nozzleAo;
//     double fShroudNdot = p * Cantera::OneAtm * fShroudVdot / Ru / Tr;
//     double fShroudMdot= fShroudNdot * fShroudM;

//    // Mass Balance
//    double outNoShroudMdot = fCoreMdot + oxCoreMdot;
//    double outNoShroudNdot = outNoShroudMdot / eqProdsMavg;

//    double outMdot = outNoShroudMdot + fShroudMdot + oxShroudMdot;
//    double outNdot = outNoShourdNdot + fShroudNdot + oxShroudNdot;

//    ndot.push_back(fCoreNdot);
//    ndot.push_back(oxCoreNdot);
//    ndot.push_back(fShroudNdot);
//    ndot.push_back(oxShroudNdot);
//    ndot.push_back(outNoShroudNdot);
//    ndot.push_back(outNdot);

//    mdot.push_back(fCoreMdot);
//    mdot.push_back(oxCoreMdot);
//    mdot.push_back(fShroudMdot);
//    mdot.push_back(oxShroudMdot);
//    mdot.push_back(outNoShroudMdot);
//    mdot.push_back(outMdot);

//    return ErrorCode::SUCCESS;
// }

// double Cfb::get_Texhaust_from_energy_balance(const std::vector<utils::KeyAnyValue>& inputs, const
// std::vector<utils::KeyAnyValue>& reacts, const std::vector<utils::KeyAnyValue>& prods, const std::vector<double>&
// nTot, double mtot, double qloss)
// {
//     std::string rxn = utils::get_value<std::string>(inputs, "reaction_mechanism");
//     double p = utils::get_value<double>(inputs, "pressure");
//     double Tr = utils::get_value<double>(inputs, "reactants_temperature");
//     double R = Cantera::GasConstant/1000.0;
//     // std::vector<double> h_prods(prods.size());
//     std::vector<double> h_reacts(reacts.size());

//     double nTotProds = nTot[2] + nTot[3] + nTot[4];
//     double nTotReacts = nTot[0] + nTot[1] + nTot[2] + nTot[3];

//     std::ostringstream rComp;
//     std::ostringstream pComp;

//     for (size_t i = 0; i < prods.size(); ++i)
//     {
//         pComp << prods[i].key << ": "  << std::any_cast<double>(prods[i].value) * nTotProds;
//         if (i == prods.size() - 1)
//             continue;
//         else
//             pComp << ", ";
//     }

//     auto sol = Cantera::newSolution(rxn);
//     auto gas = sol->thermo();

//     // Products
//     // gas->setState_TPX(Tguess, p*Cantera::OneAtm, pComp.str());

//     // gas->getEnthalpy_RT(h_prods.data());

//     // for (size_t i = 0; i < h_prods.size(); ++i)
//     // {
//         // h_prods[i] = h_prods[i] * R * Tguess;
//     // }

//     // Reactants
//     for (size_t i = 0; i < reacts.size(); ++i)
//     {
//         rComp << reacts[i].key << ": "  << std::any_cast<double>(reacts[i].value) * nTotReacts;
//         if (i == reacts.size() - 1)
//             continue;
//         else
//             rComp << ", ";
//     }
//     gas->setState_TPX(Tr, p*Cantera::OneAtm, rComp.str());
//     std::vector<double> h_RT_reacts(gas->nSpecies());
//     // std::vector<double> h_reacts(gas->nSpecies());
//     gas->getEnthalpy_RT(h_RT_reacts.data());

//     size_t speciesIdx;
//     for (size_t i = 0; i < h_reacts.size(); ++i)
//     {
//         speciesIdx = gas->speciesIndex(reacts[i].key);
//         h_reacts[i] = h_RT_reacts[speciesIdx] * R * Tr;
//     }

//     // double Eprods = 0;
//     // for (size_t i = 0; i < h_prods.size(); ++i)
//     // {
//         // Eprods += nTotProds * std::any_cast<double>(prods[i].value) * h_prods[i];
//     // }

//     double Ereacts = 0;
//     for (size_t i = 0; i < h_reacts.size(); ++i)
//     {
//         Ereacts += nTotReacts * std::any_cast<double>(reacts[i].value) * h_reacts[i];
//     }

//     double Eprods = Ereacts + qloss;

//     double hprods = Eprods/mtot;

//     gas->setState_TPX(300, p*Cantera::OneAtm, pComp.str());
//     gas->setState_HP(hprods, p);

//     double Texhaust = gas->temperature();
//     return Texhaust;
// }

} // namespace model
