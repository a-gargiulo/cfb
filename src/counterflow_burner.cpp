#include "model.h"
#include <sstream>
#include <vector>

namespace model
{

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

} // namespace model
