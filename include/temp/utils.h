#ifndef UTILS_H
#define UTILS_H

// #include "error.h"
// #include "logging.h"
// #include <any>
// #include <string>
// #include <sstream>
// #include <typeinfo>
// #include <utility>
// #include <vector>
#include <iostream>
#include <variant>

namespace utils
{

// struct KeyValue
// {
//     std::string key;
//     std::variant<int, double, std::string> value;
// };

// template <typename T>
// std::pair<T, ErrorCode> get_value(const std::vector<KeyAnyValue>& vec, const std::string& key)
// {
//     std::string errMsg;
//     logging::Logger& logger = logging::Logger::get_instance();
//     for (const KeyAnyValue& kv : vec)
//     {
//         if (kv.key == key)
//         {
//             if (kv.value.type() == typeid(T))
//                 return {std::any_cast<T>(kv.value), ErrorCode::SUCCESS};
//             else
//             {
//                 errMsg = (std::ostringstream() << "Wrong type specified for " << key << ".").str();
//                 logger.log(logging::LogLevel::ERROR, errMsg);
//                 return {T(), ErrorCode::TYPE_MISMATCH};
//             }
//         }
//     }
//     errMsg = (std::ostringstream() << key << "not found.").str();
//     logger.log(logging::LogLevel::ERROR, errMsg);
//     return {T(), ErrorCode::KEY_NOT_FOUND};
// }

} // namespace utils

#endif // UTILS_H
