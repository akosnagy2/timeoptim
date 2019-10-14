#pragma once
#include "Profile.h"

extern "C" {
#include "MTSOS\MTSOS.h"
#include "MTSOS_LP\MTSOS_LP.h"
}

class ProfileMTSOS : public Profile
{
public:
	ProfileMTSOS(const Path &path, const Robot &robot, bool linearObjective) : _linearObj(linearObjective), Profile(path, robot, "MTSOS") {};
private:
	virtual void addVelocityConstraints();
	virtual void addAccelerationConstraints();
	virtual void addDynamicConstraints();
	virtual void addObjectiveFunction();
	virtual void createModelImpl();
	virtual bool solveImpl();
	virtual void calculateTrajectoryImpl();
	bool _linearObj;
	std::vector<double> bMax, aMax;
	std::vector<double> vars;
	std::vector<double> _pathRaw;
	problem_params p_params;
	algorithm_flags a_flags;
	algorithm_params a_params;
	optional_params o_params;
};