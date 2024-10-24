#include "application.h"
#include "parsing.h"
// #include "utils.h"
// #include "model.h"
#include "logging.h"
#include <iomanip>
#include <iostream>
#include <string>

Application::Application() : m_appName("cfb"), m_isRunning(false) {};

Application::Application(const std::string& appName) : m_appName(appName), m_isRunning(false) {};

Application::~Application()
{
    shutdown();
}

void Application::print_title()
{
    const std::string version = "1.0.1";
    const std::string author = "Aldo Gargiulo";
    const std::string date = "10/24/24";

    const int headerWidth = 60;

    std::cout << std::string(headerWidth, '*') << std::endl;
    std::cout << "* " << std::setw(headerWidth - 3) << std::left << "Welcome to the '" + m_appName + "' Application"
              << "*" << std::endl;
    std::cout << "* Version " << std::setw(headerWidth - 11) << std::left << version << "*" << std::endl;
    std::cout << "* " << std::setw(headerWidth - 3) << std::left << author + ", " + date << "*" << std::endl;
    std::cout << std::string(headerWidth, '*') << std::endl;
    std::cout << std::endl;
}

void Application::initialize()
{
    Logger& logger = Logger::get_instance(LogLevel::DEBUG);
    print_title();
    m_isRunning = true;
}

bool Application::run()
{
    Logger& logger = Logger::get_instance();
    parsing::InputParser parser("input.txt");
    parsing::InputData inputs;
    inputs.reserve(19);

    if (!m_isRunning)
    {
        logger.log(LogLevel::INFO, "The application could not be initialized.");
        return true;
    }
    logger.log(LogLevel::INFO, "The application was initialized successfully.");

    if (!parser.parse(inputs))
        return true;

    auto printValue = [](const std::variant<int, double, bool, std::string>& value) {
        std::visit([](auto&& arg) {
            std::cout << arg; // Print the actual value
        }, value);
    };

    for (const auto& pair : inputs)
    {
        std::cout << "Key: " << pair.first << ", Value: ";
        printValue(pair.second); 
        std::cout << std::endl;
    }


   
    // model::Cfb cfb;
    // double Texhaust = cfb.compute_exhaust_temperature(inputs);
    // std::cout << "SOLUTION: " << Texhaust << std::endl;

    return false;
}

void Application::shutdown()
{
    Logger& logger = Logger::get_instance();
    logger.log(LogLevel::INFO, "Application is shutting down.");
}
