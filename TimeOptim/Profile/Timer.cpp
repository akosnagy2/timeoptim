#include "Timer.h"
#include <exception>
#include <windows.h>
#include <iostream>

using namespace std;

Timer::Timer()
{
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		throw exception("QueryPerformanceFrequency failed.");

	_PCfreq = double(li.QuadPart);

	QueryPerformanceCounter(&li);
	_counterStart = li.QuadPart;
}

double Timer::GetTime() const
{
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - _counterStart) / _PCfreq;
}