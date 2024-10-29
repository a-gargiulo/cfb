#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "parsing.h"
#include "logging.h"
// #include "utils.h"
#include "strutil.h"

namespace parsing
{

InputParser::InputParser(const std::string& fileName) : m_fileName(fileName){};

void InputParser::parse(InputData& data, int& err) const
{
    Logger& logger = Logger::get_instance();

    std::ifstream file(m_fileName);

    if (!file.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << m_fileName << " file.";
        logger.log(LogLevel::ERROR, oss.str());
        err = -1;
        return;
    }

    read_input_file(data, file, err);
    if (err)
    {
        std::ostringstream oss;
        oss << "Error while reading " << m_fileName << ".";
        logger.log(LogLevel::ERROR, oss.str());
        return;
    }

    file.close();

    logger.log(LogLevel::INFO, "Input file was successfully read.");

    err = 0;
    return; 
}

void InputParser::read_input_file(InputData& data, std::ifstream& file, int& err) const 
{
    std::string line;
    std::string key;
    std::string value;

    while (std::getline(file, line))
    {
        std::string_view sv = line;
        std::string_view sv_trim = strutil::trim(sv);
        if (!strutil::is_comment_or_blank(sv_trim))
        {
            auto sep = sv_trim.find('=');

            key = strutil::transform_key_string(strutil::trim(sv_trim.substr(0, sep)));
            value = strutil::trim(sv_trim.substr(sep + 1, sv_trim.size() - sep));

            if (strutil::is_bool(value))
                data[key] = strutil::to_bool(value);
            else if (strutil::is_double(value))
                data[key] = std::stod(value);
            else
                data[key] = value;
        }
    }

    err = 0;
    return;
}


}
