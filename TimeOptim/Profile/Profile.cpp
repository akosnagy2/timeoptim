#include "Profile.h"
#include <vector>
#include <numeric>
#include <iostream>

using namespace std;

bool Profile::solve()
{
	bool ok = false;

	cout << "Model solving started." << endl;
	double StartTime = _timer.GetTime();

	ok = solveImpl();

	_solveTimes.push_back(_timer.GetTime() - StartTime);

	if (ok)
		calculateTrajectory();
	return ok;
}

bool Profile::createModel(bool useDynamics)
{
	cout << "Model creation started." << endl;
	double StartTime = _timer.GetTime();

	//if (_path.getTheta() == 0) //TODO: ellenorzesek: theta, palya meret

	addVelocityConstraints();
	addAccelerationConstraints();
	if (useDynamics)
		addDynamicConstraints();

	addObjectiveFunction();

	createModelImpl();

	_modelTimes.push_back(_timer.GetTime() - StartTime);
	return true;
}

void Profile::calculateTrajectory()
{
	//Path calculation
	size_t N = _path.getLength();
	_traj.SetPath(_path);

	//Time calculation
	vector<double> time;
	time.push_back(0.0);
	for (size_t i = 1; i < N; i++)
	{
		if (abs(_a[i - 1]) < 1e-5)
			time.push_back(time.back() + _path.getTheta() / sqrt(_b[i]));
		else
			time.push_back(time.back() + (sqrt(_b[i]) - sqrt(_b[i - 1])) / _a[i - 1]);
	}
	_traj.SetTime(time);

	//Velocity, Acceleration
	calculateTrajectoryImpl();	
}

double Profile::getAvgSolveTime() const
{
	return std::accumulate(_solveTimes.begin(), _solveTimes.end(), 0.0) / (double)_solveTimes.size();
}

double Profile::getAvgModelCreationTime() const
{
	return std::accumulate(_modelTimes.begin(), _modelTimes.end(), 0.0) / (double)_modelTimes.size();
}