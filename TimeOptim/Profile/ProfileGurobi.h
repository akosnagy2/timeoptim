#pragma once
#include "Profile.h"
#include "gurobi_c++.h"

#define LP_NEW

class ProfileGurobi : public Profile
{
public:
	ProfileGurobi(const Path &path, const Robot &robot, bool linearObjective) : 
		_linearObj(linearObjective), 
		Profile(path, robot, (std::string)"Gurobi-" + (linearObjective ? "LP" : "SOCP")), _model(_env) {};
	void addFrictionConstraint(const Path &pose);
private:
	virtual void addVelocityConstraints();
	virtual void addAccelerationConstraints();
	virtual void addDynamicConstraints();
	virtual void addObjectiveFunction();
	virtual void createModelImpl();
	virtual bool solveImpl();
	virtual void calculateTrajectoryImpl();
	bool _linearObj;
	GRBEnv _env;
	GRBModel _model;
	std::vector<GRBVar> b, a, u;
	std::vector<std::vector<GRBVar>> tau;
};