#pragma once

class Timer
{
public:
	Timer();
	double GetTime() const;
private:
	double _PCfreq;
	__int64 _counterStart;
};