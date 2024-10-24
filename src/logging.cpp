#include "logging.h"

#include <iostream>

Logger& Logger::get_instance(LogLevel level)
{
    static Logger instance(level);
    return instance;
}

Logger::Logger(LogLevel level) : m_logLevel(level)
{
}

Logger::~Logger()
{
    close_log_file();
}

void Logger::set_log_level(LogLevel level)
{
    m_logLevel = level;
}

void Logger::close_log_file()
{
    if (m_logFile.is_open())
    {
        m_logFile.close();
    }
}

void Logger::log(LogLevel level, const std::string& message)
{
    if (level >= m_logLevel)
    {
        std::ostream& out = (m_logFile.is_open()) ? m_logFile : std::cout;
        out << "[" << to_string(level) << "] " << message << std::endl;
    }
}

void Logger::set_log_file(const std::string& fileName)
{
    if (m_logFile.is_open())
    {
        m_logFile.close();
    }
    m_logFile.open(fileName, std::ios::app);
}

std::string Logger::to_string(LogLevel level) const
{
    switch (level)
    {
    case LogLevel::DEBUG:
        return "DEBUG";
    case LogLevel::INFO:
        return "INFO";
    case LogLevel::WARN:
        return "WARN";
    case LogLevel::ERROR:
        return "ERROR";
    default:
        return "UNKNOWN";
    }
}
