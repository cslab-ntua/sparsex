#include "Logger.hpp"
#include "LoggerUtil.hpp"
#include <cstdio>
#include <iostream>
#include <string>

using namespace logging;

int main(int argc, char **argv)
{
    DisableLevel(Error);
    DisableLevel(Info);
    for (size_t i = 0; i < 1000000; i++) {
        LOG_ERROR << "error\n";
        LOG_INFO << "info\n";
    }

    return 0;
}
