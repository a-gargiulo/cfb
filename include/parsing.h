#ifndef PARSING_H
#define PARSING_H

#include <string>
#include <unordered_map>
#include <variant>

namespace parsing
{

using InputDataValue = std::variant<double, bool, std::string>;
using InputData = std::unordered_map<std::string, InputDataValue>;

class InputParser
{
  public:
    InputParser(const std::string& fileName);
    void parse(InputData& data, int& err) const;

  private:

    void read_input_file(InputData& data, std::ifstream& file, int& err) const;

  private:
    const std::string m_fileName;
};

} // namespace parsing

#endif // PARSING_H
