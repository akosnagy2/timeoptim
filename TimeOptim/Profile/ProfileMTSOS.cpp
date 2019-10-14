#include "ProfileMTSOS.h"
#include <numeric>
#include <iostream>
#include <algorithm>

void ProfileMTSOS::addVelocityConstraints()
{
	Path s_prime = _path.ApproximateMiddleFirstDerivative();
	double lastS = sqrt((s_prime.back() * s_prime.back()).getMaximum());

	bMax = _robot.getBConstraintForPath(_path);

	//Convert last b constraint to acceleration cosntraint	
	bMax.pop_back();
	double bMaxLast = 2 * _robot.getAccelerationConstraint()[0] * _path.getTheta() / lastS;
	if (bMaxLast < bMax.back())
	{
		bMax.pop_back();
		bMax.push_back(bMaxLast);
	}

	//Duplicate last velocity constraint
	bMax.push_back(bMax.back());
}
void ProfileMTSOS::addAccelerationConstraints()
{
	aMax = _robot.getAConstraintForPath(_path);
}
void ProfileMTSOS::addDynamicConstraints()
{

}
void ProfileMTSOS::addObjectiveFunction()
{

}

bool ProfileMTSOS::solveImpl()
{
	int N = _path.getLength();
	double* b, *u, *v, *timers;
	double fx;
	int iterations = 0;
	b = NULL;
	v = NULL;
	u = NULL;
	timers = NULL;
	fx = 0;

	if (_linearObj)
		so_MTSOS_LP(&p_params, &a_flags, &a_params, &o_params, &b, &u, &v, &fx, &timers, &iterations);
	else
		so_MTSOS(&p_params, &a_flags, &a_params, &o_params, &b, &u, &v, &fx, &timers, &iterations);

	for (int i = 0; i < N; ++i)
		_b.push_back(b[i]);
	_b[N - 1] = 0;
	for (int i = 0; i < (N - 1); ++i)
		_a.push_back((_b[i + 1] - _b[i]) / (2 * _path.getTheta()));

	if (b != NULL) {
		free(b);
	}
	if (u != NULL) {
		free(u);
	}
	if (v != NULL) {
		free(v);
	}
	if (timers != NULL) {
		free(timers);
	}

	return true;
}

void ProfileMTSOS::createModelImpl()
{
	int N = _path.getLength();
	_pathRaw = _path.getVector();

	p_params.S = _pathRaw.data();
	p_params.State_size = _path.getDimension();
	p_params.U_size = _path.getDimension();
	p_params.S_length = N;
	p_params.initial_velocity = 0;

	a_flags.timer = 1;
	a_flags.display = 0;
	a_flags.kappa = 1;

	a_params.kappa = 0;
	a_params.alpha = 0;
	a_params.beta = 0;
	a_params.epsilon = 0.00001;
	a_params.MAX_ITER = 0;

	auto pars = _robot.getParams();
	std::vector<double> torg = _robot.getTorgueConstraint();
	vars.insert(vars.end(), pars.begin(), pars.end());
	vars.insert(vars.end(), torg.begin(), torg.end());
	vars.insert(vars.end(), aMax.begin(), aMax.end());
	vars.insert(vars.end(), bMax.begin(), bMax.end());

	o_params.variables = vars.data();
	o_params.variables_length = vars.size();
	o_params.initial_b = 0;
	o_params.initial_u = 0;
}

void ProfileMTSOS::calculateTrajectoryImpl()
{
	//Acceleration calculation
	size_t N = _path.getLength();
	Path accel(_path.getDimension());
	Path s_m_prime = _path.ApproximateMiddleFirstDerivative();
	Path s_m_dprime = _path.ApproximateMiddleSecondDerivative();

	for (size_t i = 0; i < N - 1; ++i)
	{
		accel.AddPoint(s_m_prime[i] * _a[i] + s_m_dprime[i] * (_b[i] + _b[i + 1]) * 0.5);
	}
	accel.AddPoint(Point(_path.getDimension(), 0.0));	
	_traj.SetAcceleration(accel);

	//Velocity calculation
	Path s_prime = _path.ApproximateFirstDerivative();
	Path vel(_path.getDimension());
	vel.AddPoint(std::vector<double>(_path.getDimension(), 0.0));
	for (size_t i = 1; i < N; i++)
		vel.AddPoint(s_prime[i] * sqrt(_b[i]));
	_traj.SetVelocity(vel);

	//TODO: dinamika szamitasa tobbfelekeppen
	/*Path s_dprime = _path.ApproximateSecondDerivative();
	_tau.clear();
	for (int i = 1; i < N - 1; ++i)
	{
	Point m = _robot.GetDynM(_path[i], s_prime[i]);
	Point c = _robot.GetDynC(_path[i], s_prime[i], s_dprime[i]);
	Point d = _robot.GetDynD(_path[i]);
	Point t(false);
	for (int j = 0; j < _path.getDimension(); ++j)
	{
	if (m[j] * c[j] >= 0)
	t.addCoord(m[j] * _a[i - 1] + c[j] * _b[i] + d[j]);
	else
	t.addCoord(m[j] * _a[i] + c[j] * _b[i] + d[j]);
	}
	Point t(false);
	_tau.AddPoint(t);
	}*/
}