/* -*- C++ -*-
 *
 * LogPolicy.h --  Logging policy interface.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "LogPolicy.h"

void FileLogPolicy::open(const char *filename)
{
    try {
        os_->open(filename, ios_base::binary);
        if (!os_->is_open()) {
            throw ios_base::failure("");
        }
    } catch (ios_base::failure &e) {
        throw ios_base::failure("[LOGGER]: failed to open output stream, "
                                "so logging will not take place.");
    }
}

void FileLogPolicy::close()
{
    if (os_) {
        os_->close();
    }
}

void FileLogPolicy::write(const string &msg)
{
    (*os_) << msg << endl;
}
