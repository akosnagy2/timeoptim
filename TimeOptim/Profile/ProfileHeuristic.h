#pragma once
#include "Profile.h"

/*
	Represent (l <= a^T*x <= u) vector constraint 
	_coeff: a vector non zero values
	_indices: indices of non zero values
	_upperLimit: u
	_lowerLimit: l
*/
class HeuristicConstraint
{
public:
	/* Single */
	HeuristicConstraint(int index, double limit1, double limit2) :
		_coeff(std::vector<double>(1, 1.0)),
		_indices(std::vector<int>(1, index)),
		_upperLimit((limit1 > limit2) ? limit1 : limit2),
		_lowerLimit((limit1 > limit2) ? limit2 : limit1) {}
	/* General */
	HeuristicConstraint(std::vector<double> coeff, std::vector<int> indices, double limit1, double limit2) :
		_coeff(coeff),
		_indices(indices),
		_upperLimit((limit1 > limit2) ? limit1 : limit2),
		_lowerLimit((limit1 > limit2) ? limit2 : limit1) {}

	std::vector<double> _coeff;
	std::vector<int> _indices;
	double _upperLimit;
	double _lowerLimit;
};

class ProfileHeuristic : public Profile
{
public:
	ProfileHeuristic(const Path &path,const Robot &robot) : Profile(path, robot, "Heuristic") {};
private:
	virtual void addVelocityConstraints();
	virtual void addAccelerationConstraints();
	virtual void addDynamicConstraints();
	virtual void addObjectiveFunction();
	HeuristicConstraint generateDynamicConstraint(int i, double m, double c, double limit1, double limit2);
	virtual void createModelImpl();
	virtual bool solveImpl();
	virtual void calculateTrajectoryImpl();
	void reduceStartStop();
	void sortConstraints();
	void heuristic4();	
	bool test(double b0, double b1, double b2, int i);
	std::vector<double> _bb;
	std::vector<HeuristicConstraint> _velConsts;
	std::vector<HeuristicConstraint> _accelConsts;
	std::vector<HeuristicConstraint> _dynConsts;
	std::vector<std::vector<HeuristicConstraint>> _sortedConsts;
};