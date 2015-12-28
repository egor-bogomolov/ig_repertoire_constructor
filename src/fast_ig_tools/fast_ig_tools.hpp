#pragma once

#include "standard_include.hpp"


std::string join_cmd_line(size_t argc, char **argv) {
    std::string result = argv[0];
    for (size_t i = 1; i < argc; ++i) {
        result += " ";
        result += argv[i];
    }

    return result;
}


void create_console_logger(std::string log_props_file = "") {
    using namespace logging;

    logger *lg = create_logger(log_props_file);
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

// vim: ts=4:sw=4
