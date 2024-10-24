#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "error.h"
#include "utils.h"
#include "model.h"
#include <unordered_map>
#include <vector>
#include <utility>

namespace model
{
class Cfb
{
  public:
    virtual model::ThermoState& compute_exhaust_temperature(const std::vector<utils::KeyAnyValue>& inputs);

  private:
    double get_exhaust_density(const std::vector<utils::KeyAnyValue>& inputs, const std::vector<utils::KeyAnyValue>& prods,
                               double& Texhaust, const std::vector<double>& nTot);

    double get_qloss_from_heat_transfer(const std::vector<utils::KeyAnyValue>& inputs,
                                        const std::vector<utils::KeyAnyValue>& reacts,
                                        const std::vector<utils::KeyAnyValue>& prods, const std::vector<double>& nTot,
                                        double Tguess);

    bool fuel_is_hydrocarbon(const std::string& fuel, std::unordered_map<std::string, int>& elementCounts);

    ErrorCode determine_reactants_composition(const std::vector<utils::KeyAnyValue>& inputs, std::unordered_map<std::string, double>& reacts, std::string& reactsStr);

     ErrorCode compute_eq_products_composition(const std::vector<utils::KeyAnyValue>& inputs,
                                                        const std::string& reactsStr,
        std::unordered_map<std::string, double>& prods                                              );

    double get_Texhaust_from_energy_balance(const std::vector<utils::KeyAnyValue>& inputs,
                                            const std::vector<utils::KeyAnyValue>& reacts,
                                            const std::vector<utils::KeyAnyValue>& prods, const std::vector<double>& nTot,
                                            double mtot, double qloss);

    ErrorCode flow_rates_from_mass_balance(const std::vector<utils::KeyAnyValue>& inputs,
                                           const std::vector<utils::KeyAnyValue>& prods,
                                           std::vector<double>& ndot,
                                           std::vector<double>& mdot);

    ErrorCode adjust_products_mole_fractions(const std::vector<utils::KeyAnyValue>& inputs,
                                             const std::vector<double>& ndot,
                                             const std::vector<utils::KeyAnyValue>& prods,
                                             std::vector<utils::KeyAnyValue>& prods_adj);

    ErrorCode adjust_reactants_mole_fractions(const std::vector<utils::KeyAnyValue>& inputs,
                                              const std::vector<double>& ndot,
                                              const std::string& rComp,
                                              std::vector<utils::KeyAnyValue>& reacts_adj);
};

} // namespace model

#endif // MODEL_H
