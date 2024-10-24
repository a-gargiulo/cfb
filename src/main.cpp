#include "application.h"

int main(int argc, char* argv[])
{
    Application app("Counterflow Burner");

    app.initialize();

    app.run();

    return 0;
}
