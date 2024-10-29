#ifndef APPLICATION_H
#define APPLICATION_H

#include <string>

class Application
{
  public:
    Application();

    Application(const std::string& name);

    ~Application();

    void initialize();

    bool run() const;

    void shutdown() const;

  private:
    void print_title() const;

  private:
    const std::string m_appName;
    bool m_isRunning;
};

#endif // APPLICATION_H
