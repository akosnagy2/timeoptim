#pragma once
#include "Path\Path.h"
#include "Path\Trajectory.h"
#include "Robot\Robot.h"
#include "Profile\Timer.h"

class Profile
{
public:
	Profile(const Path &path, const Robot &robot, std::string name) :
		_traj(path.getDimension()),
		_path(path),
		_robot(robot),
		_name(name) {};
	const std::string& getName() const { return _name; }
	bool solve();
	bool createModel(bool useDynamics);
	const std::vector<double>& getSolveTimes() const { return _solveTimes; }
	const std::vector<double>& getModelCreationTimes() const { return _modelTimes; }
	double getAvgSolveTime() const;
	double getAvgModelCreationTime() const;
	const Trajectory& getTrajectory() const { return _traj; }
	const std::vector<double>& getB() const { return _b; }
protected:
	void calculateTrajectory();
	virtual void calculateTrajectoryImpl() = 0;
	virtual void addVelocityConstraints() = 0;
	virtual void addAccelerationConstraints() = 0;
	virtual void addDynamicConstraints() = 0;
	virtual void addObjectiveFunction() = 0;
	virtual void createModelImpl() = 0;
	virtual bool solveImpl() = 0;
	//Info
	Timer _timer;
	std::string _name;
	std::vector<double> _solveTimes, _modelTimes;
	//Inputs
	const Robot& _robot;
	const Path& _path;
	//Outputs
	Trajectory _traj;
	std::vector<double> _b, _a;
};