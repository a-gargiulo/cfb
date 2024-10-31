#include "application.h"
#include "parsing.h"
// #include "utils.h"
#include "model.h"
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

void Application::print_title() const
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
    logger.log(LogLevel::DEBUG, "Logger succesfuly initialized.");
    print_title();
    m_isRunning = true;
}

bool Application::run() const
{
    int err;
    Logger& logger = Logger::get_instance();

    parsing::InputParser parser("input.txt");
    parsing::InputData inputs;
    inputs.reserve(20);

    if (!m_isRunning)
    {
        logger.log(LogLevel::INFO, "The application could not be initialized.");
        return true;
    }
    logger.log(LogLevel::INFO, "The application was initialized successfully.");

    parser.parse(inputs, err);
    if (err)
        return true;

   
    model::Equilibrium cfb(inputs);
    cfb.compute_exhaust_thermo_state(inputs, err);
    // std::cout << "SOLUTION: " << Texhaust << std::endl;

    return false;
}

void Application::shutdown() const
{
    Logger& logger = Logger::get_instance();
    logger.log(LogLevel::INFO, "Application is shutting down.");
}


    // auto printValue = [](const std::variant<int, double, bool, std::string>& value) {
    //     std::visit([](auto&& arg) {
    //         std::cout << arg;
    //     }, value);
    // };

    // for (const auto& pair : inputs)
    // {
    //     std::cout << "Key: " << pair.first << ", Value: ";
    //     printValue(pair.second); 
    //     std::cout << std::endl;
    // }
