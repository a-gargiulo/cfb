#include <algorithm>
#include <cctype>
#include <ios>
#include <sstream>
#include <string>
#include <string_view>
#include <map>

#include "strutil.h"

std::string_view strutil::trim(const std::string_view sv)
{
    auto start = sv.find_first_not_of(" \t\n\v\f\r\"\'");
    if (start == std::string_view::npos)
        return {};

    auto end = sv.find_last_not_of(" \t\n\v\f\r\"\'");

    return sv.substr(start, end - start + 1);
}

bool strutil::is_comment_or_blank(const std::string_view sv)
{
    return (sv.substr(0, 2) == "//" || sv == "");
}

std::string strutil::transform_key_string(std::string_view sv)
{
    std::string_view::size_type paren = sv.find("(");

    std::string str(strutil::trim(sv.substr(0, paren)));

    std::transform(str.begin(), str.end(), str.begin(),
                   [](unsigned char c) { return std::isspace(c) ? '_' : std::tolower(c); });

    return str;
}

bool strutil::is_bool(const std::string& str)
{
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return (lower_str == "true" || lower_str == "false" || lower_str == "t" || lower_str == "f");
}

bool strutil::is_double(const std::string& str)
{
    std::istringstream iss(str);
    double d;
    iss >> std::noskipws >> d;
    return iss.eof() && !iss.fail();
}

bool strutil::to_bool(const std::string& str)
{
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return lower_str == "true" || lower_str == "t";
}

bool strutil::is_digit(char ch)
{
    return std::isdigit(static_cast<unsigned char>(ch));
}


// std::vector<std::pair<std::string, double>> strutil::process_composition_string(const std::string& str) {
//     std::vector<std::pair<std::string, double>> speciesMap;
//     std::istringstream iss(str);
//     std::string segment;

//     // Split the str string by commas
//     while (std::getline(iss, segment, ',')) {
//         std::istringstream segmentStream(segment);
//         std::string species;
//         double value;

//         // Split each segment by the colon to get the species and its value
//         if (std::getline(segmentStream, species, ':') && segmentStream >> value) {
//             // Remove any leading/trailing whitespace from the species string
//             species.erase(0, species.find_first_not_of(" \t"));
//             species.erase(species.find_last_not_of(" \t") + 1);

//             speciesMap.push_back({species, value});
//         }
//     }

//     return std::move(speciesMap);
// }
