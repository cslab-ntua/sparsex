/*
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file TimerTest.cpp
 * \brief Test timer
 *
 * \author Vasileios Karakasis
 * \date 2013
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Timer.hpp>
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace sparsex::timing;

int main()
{
    Timer timer;

	cout << "sleeping for 1 second and reporting:\n";
	timer.Start();
	sleep(1);
	timer.Pause();
	cout << "timed %lf secs\n", timer.ElapsedTime();

    timer.Clear();
	cout << "sleeping 2 times for 1 second and reporting:\n";
	timer.Start();
	sleep(1);
	timer.Pause();

	sleep(1);

	timer.Start();
	sleep(1);
	timer.Pause();
	cout << "timed %lf secs\n", timer.ElapsedTime();
	
	return 0;
}
