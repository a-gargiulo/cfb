#ifndef STRUTIL_H
#define STRUTIL_H

#include <string>
#include <string_view>
#include <vector>
#include <utility>

namespace strutil
{
std::string_view trim(const std::string_view sv);
bool is_comment_or_blank(const std::string_view sv);
std::string transform_key_string(std::string_view sv);
bool is_double(const std::string& str);
bool is_bool(const std::string& str);
bool to_bool(const std::string& str);
bool is_digit(char ch);
std::vector<std::pair<std::string, double>> process_composition_string(const std::string& str);
} // namespace strutil

#endif // STRUTIL_H
