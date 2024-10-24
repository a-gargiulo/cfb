#ifndef LOGGING_H
#define LOGGING_H

#include <fstream>
#include <string>


enum class LogLevel
{
    DEBUG,
    INFO,
    WARN,
    ERROR
};

class Logger
{
  public:
    static Logger& get_instance(LogLevel level = LogLevel::INFO);

    ~Logger();

    void set_log_level(LogLevel level);

    void log(LogLevel level, const std::string& message);
    void set_log_file(const std::string& fileName);

  private:
    Logger(LogLevel level);

    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    Logger(Logger&&) = delete;
    Logger& operator=(Logger&&) = delete;

    LogLevel m_logLevel;
    std::ofstream m_logFile;

    void close_log_file();
    std::string to_string(LogLevel level) const;
};

#endif // LOGGING_H
