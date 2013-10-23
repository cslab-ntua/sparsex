/* -*- C++ -*-
 *
 * OutputPolicy.cc --  Output policies for logging.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "OutputPolicy.hpp"

namespace logging {

void StdPolicy::operator()(const buffer_type &buffer)
{
    if (pthread_mutex_lock(&write_mutex_) != 0) {
        cerr << "[LOGGER]: failed to lock shared variable." 
             << endl;
        exit(1);
    }

    cout << buffer.str();

    if (pthread_mutex_unlock(&write_mutex_) != 0) {
        cerr << "[LOGGER]: failed to unlock shared variable." 
             << endl;
        exit(1);
    }
}

void FilePolicy::Init(const char *log_file)
{
    try {
        os_.open(log_file, ios_base::binary);
        if (!os_.is_open()) {
            throw ios_base::failure("");
        }
    } catch (ios_base::failure &e) {
        throw ios_base::failure("[LOGGER]: failed to open output stream, "
                                "so logging will not take place.");
    }
}

void FilePolicy::operator()(const buffer_type &buffer)
{
    if (pthread_mutex_lock(&write_mutex_) != 0) {
        cerr << "[LOGGER]: failed to lock shared variable." 
             << endl;
        exit(1);
    }

    os_ << buffer.str();

    if (pthread_mutex_unlock(&write_mutex_) != 0) {
        cerr << "[LOGGER]: failed to unlock shared variable." 
             << endl;
        exit(1);
    }
}

}
