#include "ProfileGurobi.h"
#include <numeric>
#include <chrono>
#include <iostream>
#include <algorithm>

using namespace std;

void ProfileGurobi::createModelImpl()
{

}

void ProfileGurobi::addVelocityConstraints()
{
	int N = (int)_path.getLength();

	/* Velocity limits */
	std::vector<double> velMax = _robot.getVelocityConstraint();

	/* Velocity variables */
	for (int i = 0; i < N; ++i)
	{
		b.push_back(_model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
	}
	_model.update();

	/* Velocity constraints */
	auto s_prime = _path.ApproximateFirstDerivative();
	for (int i = 0; i < N; ++i)
	{
		if ((i == 0) || (i == (N - 1)))
			_model.addConstr(b[i] == 0);
		else
		{
			for (int j = 0; j < _path.getDimension(); ++j)
			{
				_model.addConstr(((s_prime[i][j] * s_prime[i][j]) * b[i]) <= (velMax[j] * velMax[j]));
			}
		}
	}
	_model.update();
}
void ProfileGurobi::addAccelerationConstraints()
{
	int N = (int)_path.getLength();

	/* Acceleration variables */
	for (int i = 0; i < N - 1; ++i)
	{
		a.push_back(_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
	}
	_model.update();

	/* Acceleration-velocity constraints */
	for (int i = 0; i < N - 1; ++i)
	{
		_model.addConstr(2 * a[i] * _path.getTheta() == (b[i + 1] - b[i]));
	}
	_model.update();

	/* Acceleration limits */
	Point accelMax = _robot.getAccelerationConstraint();

	/* Acceleration constraints */
	auto s_prime = _path.ApproximateFirstDerivative();
	auto s_dprime = _path.ApproximateSecondDerivative();
	for (int i = 1; i < N - 1; ++i)
	{
		for (int j = 0; j < _path.getDimension(); ++j)
		{
			if (s_prime[i][j] * s_dprime[i][j] >= 0)
			{
				_model.addConstr((s_prime[i][j] * a[i-1] + s_dprime[i][j] * b[i]) <= accelMax[j]);
				_model.addConstr((s_prime[i][j] * a[i-1] + s_dprime[i][j] * b[i]) >= -accelMax[j]);
			}
			else
			{
				_model.addConstr((s_prime[i][j] * a[i] + s_dprime[i][j] * b[i]) <= accelMax[j]);
				_model.addConstr((s_prime[i][j] * a[i] + s_dprime[i][j] * b[i]) >= -accelMax[j]);
			}
		}
	}
	_model.update();

}
void ProfileGurobi::addDynamicConstraints()
{
	size_t N = _path.getLength();

	/* Dynamics constraints */
	for (int i = 0; i < N - 2; ++i)
	{
		vector<GRBVar> e;
		for (int j = 0; j < _path.getDimension(); ++j)
		{
			e.push_back(_model.addVar(-_robot.getUConstraintForPath(_path)[i][j], _robot.getUConstraintForPath(_path)[i][j], 0.0, GRB_CONTINUOUS));
		}
		tau.push_back(e);
	}
	_model.update();
	Path fd = _path.ApproximateFirstDerivative();
	Path sd = _path.ApproximateSecondDerivative();
	for (int i = 1; i < N - 1; ++i)
	{
		Point m = _robot.GetDynM(_path[i], fd[i]);
		Point c = _robot.GetDynC(_path[i], fd[i], sd[i]);
		Point d = _robot.GetDynD(_path[i]);
		for (int j = 0; j < _path.getDimension(); ++j)
		{
			if (m[j] * c[j] >= 0)
				_model.addConstr(m[j] * a[i - 1] + c[j] * b[i] + d[j] == tau[i - 1][j]);
			else
				_model.addConstr(m[j] * a[i] + c[j] * b[i] + d[j] == tau[i - 1][j]);
		}
	}
	_model.update();
}

void ProfileGurobi::addFrictionConstraint(const Path &pose)
{
	size_t N = pose.getLength();

	/* Acceleration limits */
	double a_max = 0.5;

	/* Acceleration constraints */
	auto r_m_prime = pose.ApproximateMiddleFirstDerivative();
	auto r_m_dprime = pose.ApproximateMiddleSecondDerivative();
	for (int i = 0; i < N - 1; ++i)
	{
		auto x = r_m_prime[i][0] * a[i] + r_m_dprime[i][0] * (b[i] + b[i + 1]) * 0.5;
		auto y = r_m_prime[i][1] * a[i] + r_m_dprime[i][1] * (b[i] + b[i + 1]) * 0.5;
		auto z = r_m_prime[i][2] * a[i] + r_m_dprime[i][2] * (b[i] + b[i + 1]) * 0.5;
		_model.addQConstr(x*x + y*y + z*z <= a_max*a_max);
	}
	_model.update();
}

void ProfileGurobi::addObjectiveFunction()
{
	size_t N = _path.getLength();

	/*  Objective function */
	GRBLinExpr objExpr;
	if (_linearObj)
	{
		for (size_t i = 0; i < b.size(); ++i)
			objExpr += b[i];

		_model.setObjective(objExpr, GRB_MAXIMIZE);
		_model.getEnv().set(GRB_IntParam_Method, 1);
		_model.getEnv().set(GRB_DoubleParam_BarQCPConvTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_FeasibilityTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_IntFeasTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_MIPGap, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_OptimalityTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_PSDTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_PerturbValue, 10e-08);
	}
	else
	{
		_model.getEnv().set(GRB_IntParam_Method, 1);
		std::vector<GRBVar> c, d, t1, t2, t3, t4;

		for (int i = 0; i < N - 1; ++i)
		{
			d.push_back(_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
			t1.push_back(_model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
			t2.push_back(_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
		}

		for (int i = 0; i < N; ++i)
		{
			c.push_back(_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
			t3.push_back(_model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
			t4.push_back(_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
		}
		_model.update();

		for (int i = 0; i < N - 1; ++i)
		{
			_model.addConstr(c[i + 1] + c[i] + d[i] == t1[i]);
			_model.addConstr(c[i + 1] + c[i] - d[i] == t2[i]);
			_model.addQConstr((4 + t2[i] * t2[i]) <= t1[i] * t1[i]);
		}

		for (int i = 0; i < N; ++i)
		{
			_model.addConstr((b[i] + 1) == t3[i]);
			_model.addConstr((b[i] - 1) == t4[i]);
			_model.addQConstr((4 * c[i] * c[i] + t4[i] * t4[i]) <= t3[i] * t3[i]);
		}

		for (int i = 0; i < N - 1; ++i)
			objExpr += 2 * d[i] * _path.getTheta();
		_model.setObjective(objExpr, GRB_MINIMIZE);
		_model.getEnv().set(GRB_DoubleParam_BarQCPConvTol, 10e-10);
		_model.getEnv().set(GRB_DoubleParam_FeasibilityTol, 10e-10);
		_model.getEnv().set(GRB_DoubleParam_IntFeasTol, 10e-10);
		_model.getEnv().set(GRB_DoubleParam_MIPGap, 10e-10);
		_model.getEnv().set(GRB_DoubleParam_OptimalityTol, 10e-10);
		_model.getEnv().set(GRB_DoubleParam_PSDTol, 10e-08);
		_model.getEnv().set(GRB_DoubleParam_PerturbValue, 10e-08);
	}
}

bool ProfileGurobi::solveImpl()
{
	size_t N = _path.getLength();

	/* Solve program */
	_model.getEnv().set(GRB_IntParam_OutputFlag, 0);
	_model.optimize();

	int status = _model.get(GRB_IntAttr_Status);

	if ((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) ||
		(status == GRB_UNBOUNDED))
	{
		cout << "The _model cannot be solved because it is "
			<< "infeasible or unbounded" << endl;
		return false;
	}

	if (status != GRB_OPTIMAL)
	{
		cout << "Optimization was stopped with status " << status << endl;
	}
	/*else
	{
		cout << "Optimal solution found." << endl;	
	}*/

	for (int i = 0; i < N; ++i)
		_b.push_back(b[i].get(GRB_DoubleAttr_X));
	for (int i = 0; i < (N - 1); ++i)
		_a.push_back(a[i].get(GRB_DoubleAttr_X));

	return true;
}

void ProfileGurobi::calculateTrajectoryImpl()
{
	//Acceleration calculation
	size_t N = _path.getLength();
	Path accel(_path.getDimension());
	if (_linearObj)
	{
		auto s_prime = _path.ApproximateFirstDerivative();
		auto s_dprime = _path.ApproximateSecondDerivative();
		accel.AddPoint(Point(_path.getDimension(), 0.0));
		for (int i = 1; i < N - 1; ++i)
		{
			Point t;
			for (int j = 0; j < _path.getDimension(); ++j)
			{
				if (s_prime[i][j] * s_dprime[i][j] > 0)
				{
					t.addCoord((s_prime[i][j] * _a[i - 1] + s_dprime[i][j] * _b[i]));
				}
				else
				{
					t.addCoord((s_prime[i][j] * _a[i] + s_dprime[i][j] * _b[i]));
				}
			}
			accel.AddPoint(t);
		}
		accel.AddPoint(Point(_path.getDimension(), 0.0));
	}
	else
	{
		Path s_m_prime = _path.ApproximateMiddleFirstDerivative();
		Path s_m_dprime = _path.ApproximateMiddleSecondDerivative();

		for (size_t i = 0; i < N - 1; ++i)
		{
			accel.AddPoint(s_m_prime[i] * _a[i] + s_m_dprime[i] * (_b[i] + _b[i + 1]) * 0.5);
		}
		accel.AddPoint(Point(_path.getDimension(), 0.0));
	}
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