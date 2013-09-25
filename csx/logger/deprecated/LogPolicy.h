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
#ifndef LIBCSX_LOGPOLICY_H__
#define LIBCSX_LOGPOLICY_H__

#include <iostream>
#include <fstream>

using namespace std;

/* Interface of a logging policy */
class LogPolicy
{
public:
    virtual void open(const char *filename) = 0;
    virtual void close() = 0;
    virtual void write(const string &msg) = 0;
};

/* Logging policy that directs all the output of the logger to stdout */ 
class StdLogPolicy : public LogPolicy
{
public:
    StdLogPolicy() {}
    ~StdLogPolicy() {}
    void open(const char *filename) {}
    void close() {}
    void write(const string &msg)
    {
        cout << msg << endl;
    }
};

/* Logging policy that directs all the output of the logger to a file */ 
class FileLogPolicy : public LogPolicy
{
public:

    FileLogPolicy()
    {
        os_ = new ofstream;
    }

    ~FileLogPolicy()
    {
        if (os_) {
            close();
            delete os_;
        }
    }
    void open(const char *filename);
    void close();
    void write(const string &msg);
private:
    ofstream *os_;
};

#endif  // LIBCSX_LOGPOLICY_H__
