#include "model.h"
#include "strutil.h"
#include <cantera/core.h>
#include <sstream>
#include <stack>
#include <unordered_map>
#include <vector>
#include <utility>

namespace model
{

CounterflowBurner::CounterflowBurner(const parsing::InputData& inputs)
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

void CounterflowBurner::process_composition_string(const std::string& str, std::vector<MoleFraction>& comp)
{
    std::istringstream iss(str);
    std::string segment;

    // Split the str string by commas
    while (std::getline(iss, segment, ','))
    {
        std::istringstream segmentStream(segment);

        std::string name;
        double fraction;

        // Split each segment by the colon to get the species name and its value
        if (std::getline(segmentStream, name, ':') && segmentStream >> fraction)
        {
            // Remove any leading/trailing whitespace from the species string
            name.erase(0, name.find_first_not_of(" \t"));
            name.erase(name.find_last_not_of(" \t") + 1);

            comp.emplace_back(name, fraction);
        }
    }
}

bool CounterflowBurner::fuel_is_hydrocarbon(const std::string& fuel, std::unordered_map<std::string, int>& nElements)
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

            nElements[element] += count;
            i = j;
        }
        else if (ch == '(')
        {
            // Push current element counts onto the stack and start a new scope
            parenStack.push(std::move(nElements));
            nElements.clear();
            i++;
        }
        else if (ch == ')')
        {
            // Pop the stack and handle the multiplier
            if (parenStack.empty())
            {
                return false; // Unmatched parenthesis
            }

            std::unordered_map<std::string, int> innerElements;
            innerElements.reserve(3);
            innerElements = std::move(nElements);
            nElements = std::move(parenStack.top());
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
            for (const auto& elm : innerElements)
            {
                nElements[elm.first] += elm.second * multiplier;
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

} // namespace model
