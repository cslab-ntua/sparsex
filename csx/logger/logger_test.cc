#include "Logger.hpp"
#include <cstdio>
#include <iostream>
#include <string>

using namespace logging;
//DISABLE_LOGGING_LEVEL(Error);

int main(int argc, char **argv)
{
    // log<logging::Error>() << "error\n";
    // log<logging::Info>() << "info\n";
//    LOG_ERROR << "error\n";
    LOG_INFO << "info\n";

    return 0;
}
