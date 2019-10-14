#include "ProfileHeuristic.h"
#include <numeric>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include "Simplex\Simplex.h"

using namespace std;

void ProfileHeuristic::createModelImpl()
{
	sortConstraints();
}

void ProfileHeuristic::addVelocityConstraints()
{
	size_t N = _path.getLength();

	/* Velocity limits */
	auto _c = _robot.getBConstraintForPath(_path);
	for (int i = 0; i < N; ++i)
		_velConsts.push_back(HeuristicConstraint(i, 0, _c[i]));
}

void ProfileHeuristic::addAccelerationConstraints()
{
	size_t N = _path.getLength();

	/* Acceleration limits */
	auto s_prime = _path.ApproximateFirstDerivative();
	auto s_dprime = _path.ApproximateSecondDerivative();
	for (int i = 1; i < N - 1; ++i)
	{
		Point gNumU = _robot.getAccelerationConstraint();
		Point gNumL = _robot.getAccelerationConstraint() * (-1.0);
		Point mi = s_prime[i] / (2 * _path.getTheta());
		Point ci = s_dprime[i];
		for (int j = 0; j < _path.getDimension(); ++j)
		{
			//assert(mi[j] != 0.0);
			//assert(ci[j] != 0.0);
			_accelConsts.push_back(generateDynamicConstraint(i, mi[j], ci[j], gNumL[j], gNumU[j]));
		}
	}
}

void ProfileHeuristic::addDynamicConstraints()
{
	size_t N = _path.getLength();

	/* Dynamic limits */
	Path fd = _path.ApproximateFirstDerivative();
	Path sd = _path.ApproximateSecondDerivative();
	for (int i = 1; i < N - 1; ++i)
	{
		Point gNumU = _robot.getTorgueConstraint() - _robot.GetDynD(_path[i]);
		Point gNumL = _robot.getTorgueConstraint() * (-1.0) - _robot.GetDynD(_path[i]);
		Point mi = _robot.GetDynM(_path[i], fd[i]) / (2 * _path.getTheta());
		Point ci = _robot.GetDynC(_path[i], fd[i], sd[i]);
		for (int j = 0; j < _path.getDimension(); ++j)
		{
			_dynConsts.push_back(generateDynamicConstraint(i, mi[j], ci[j], gNumL[j], gNumU[j]));
		}
	}
}

HeuristicConstraint ProfileHeuristic::generateDynamicConstraint(int i, double m, double c, double limit1, double limit2)
{
	vector<double> a;
	vector<int> ind;
	if (m * c >= 0)
	{
		a.push_back(-m);
		a.push_back(c + m);
		ind.push_back(i - 1);
		ind.push_back(i);

	}
	else
	{
		a.push_back(-(m - c));
		a.push_back(m);
		ind.push_back(i);
		ind.push_back(i + 1);
	}
	return HeuristicConstraint(a, ind, limit1, limit2);
}

void ProfileHeuristic::addObjectiveFunction()
{

}

bool ProfileHeuristic::solveImpl()
{
	size_t N = _path.getLength();

	reduceStartStop();
	heuristic4();

	for (int i = 0; i < (N - 1); ++i)
		_a.push_back((_b[i + 1] - _b[i]) / (2 * _path.getTheta()));

	return true;
}

void ProfileHeuristic::sortConstraints()
{
	_sortedConsts.resize(_path.getLength() - 1);
	for (size_t i = 0; i < _accelConsts.size(); ++i)
	{
		_sortedConsts[_accelConsts[i]._indices[0]].push_back(_accelConsts[i]);
	}
	for (size_t i = 0; i < _dynConsts.size(); ++i)
	{
		_sortedConsts[_dynConsts[i]._indices[0]].push_back(_dynConsts[i]);
	}

	for (size_t i = 0; i < _sortedConsts.size(); ++i)
	{
		auto iter = _sortedConsts[i].begin();
		while (iter != _sortedConsts[i].end())
		{
			assert(iter->_indices.size() == 2);
			if ((iter->_coeff[0] == 0.0) && (iter->_coeff[1] == 0.0))
			{
				iter = _sortedConsts[i].erase(iter);
			}
			else if (iter->_coeff[0] == 0.0)
			{
				double m;
				if (iter->_coeff[1] > 0.0)
					m = iter->_upperLimit / iter->_coeff[1];
				else
					m = iter->_lowerLimit / iter->_coeff[1];				
				_velConsts[iter->_indices[1]]._upperLimit = min(_velConsts[iter->_indices[1]]._upperLimit, m);
				iter = _sortedConsts[i].erase(iter);
			}
			else if (iter->_coeff[1] == 0.0)
			{
				double m;
				if (iter->_coeff[0] > 0.0)
					m = iter->_upperLimit / iter->_coeff[0];
				else
					m = iter->_lowerLimit / iter->_coeff[0];

				_velConsts[iter->_indices[0]]._upperLimit = min(_velConsts[iter->_indices[0]]._upperLimit, m);
				iter = _sortedConsts[i].erase(iter);
			}					
			else
			{
				++iter;
			}
		}
	}
}

void ProfileHeuristic::reduceStartStop()
{
	size_t N = _path.getLength();

	_b.resize(N);
	_b[0] = 0;
	_b[N - 1] = 0;

	/* Delete first, last limits */
	_velConsts.erase(_velConsts.begin());
	_velConsts.pop_back();

	/* Change limits*/
	for (size_t i = 0; i < _sortedConsts[0].size(); ++i)
	{
		if (_sortedConsts[0][i]._coeff[1] > 0.0)
			_velConsts[0]._upperLimit = min(_velConsts[0]._upperLimit, _sortedConsts[0][i]._upperLimit / (_sortedConsts[0][i]._coeff[1]));
		else
			_velConsts[0]._upperLimit = min(_velConsts[0]._upperLimit, _sortedConsts[0][i]._lowerLimit / (_sortedConsts[0][i]._coeff[1]));
	}
	for (size_t i = 0; i < _sortedConsts[N - 2].size(); ++i)
	{
		if (_sortedConsts[N - 2][i]._coeff[0] < 0.0)
			_velConsts[N - 3]._upperLimit = min(_velConsts[N - 3]._upperLimit, _sortedConsts[N - 2][i]._lowerLimit / (_sortedConsts[N - 2][i]._coeff[0]));
		else
			_velConsts[N - 3]._upperLimit = min(_velConsts[N - 3]._upperLimit, _sortedConsts[N - 2][i]._upperLimit / (_sortedConsts[N - 2][i]._coeff[0]));

	}
	/* Delete first, last limits */
	_sortedConsts.erase(_sortedConsts.begin());
	_sortedConsts.pop_back();
}


void ProfileHeuristic::heuristic4()
{
	Tableau tab;
	size_t N = _path.getLength() - 2;

	_bb.resize(N);

	/* Init at maximum velocity */
	for (int j = 0; j < N; ++j)
		_bb[j] = _velConsts[j]._upperLimit;	

	for (int i = 1; i < N; ++i)
	{	
		//Solve 1D LP problem
		double b_max = _bb[i];
		double b_min = 0;
		for (size_t j = 0; j < _sortedConsts[i - 1].size(); ++j)
		{
			b_max = min(b_max, (-_sortedConsts[i - 1][j]._coeff[0] * _bb[i - 1] + _sortedConsts[i - 1][j]._upperLimit)/ _sortedConsts[i - 1][j]._coeff[1]);
			b_min = max(b_min, (_sortedConsts[i - 1][j]._lowerLimit - _bb[i - 1] * _sortedConsts[i - 1][j]._coeff[0])/ _sortedConsts[i - 1][j]._coeff[1]);
		}			
		if (b_min <= b_max) //1D LP
		{
			_bb[i] = b_max;
		}		
		else //2D LP
		{
			tab.m = 3 + 2 * _sortedConsts[i - 1].size();
			tab.n = 3;

			/* Objective Fcn*/
			tab.mat[0][0] = 0;
			tab.mat[0][1] = -1;
			tab.mat[0][2] = -1;
			for (int j = 3; j < SIMPLEX_M; ++j)
				tab.mat[0][j] = 0;

			/* Constraints */
			/* Velocity */
			tab.mat[1][0] = _bb[i - 1];
			tab.mat[1][1] = 1;
			tab.mat[1][2] = 0;

			tab.mat[2][0] = _bb[i];
			tab.mat[2][1] = 0;
			tab.mat[2][2] = 1;

			/* Acceleration, dynamics */
			for (size_t j = 0; j < _sortedConsts[i - 1].size(); ++j)
			{
				assert(abs(_sortedConsts[i - 1][j]._upperLimit) > 1e-08);
				tab.mat[3 + j][0] = _sortedConsts[i - 1][j]._upperLimit;
				tab.mat[3 + j][1] = _sortedConsts[i - 1][j]._coeff[0];
				tab.mat[3 + j][2] = _sortedConsts[i - 1][j]._coeff[1];
			}
			for (size_t j = 0; j < _sortedConsts[i - 1].size(); ++j)
			{
				assert(abs(_sortedConsts[i - 1][j]._upperLimit) > 1e-08);
				tab.mat[3 + _sortedConsts[i - 1].size() + j][0] = -_sortedConsts[i - 1][j]._lowerLimit;
				tab.mat[3 + _sortedConsts[i - 1].size() + j][1] = -_sortedConsts[i - 1][j]._coeff[0];
				tab.mat[3 + _sortedConsts[i - 1].size() + j][2] = -_sortedConsts[i - 1][j]._coeff[1];
			}

			simplex(&tab);
			double bb[2];
			get_optimal_vector(&tab, bb);
			_bb[i] = bb[1];
		}
	}
	for (int i = N - 2; i >= 0; --i)
	{
		for (size_t j = 0; j < _sortedConsts[i].size(); ++j)
		{
			if (_sortedConsts[i][j]._coeff[0] < 0.0)
				_bb[i] = min(_bb[i], (_sortedConsts[i][j]._coeff[1]*_bb[i+1] - _sortedConsts[i][j]._lowerLimit)/(-_sortedConsts[i][j]._coeff[0]));
			else
				_bb[i] = min(_bb[i], (-_sortedConsts[i][j]._coeff[1] * _bb[i + 1] + _sortedConsts[i][j]._upperLimit) / (_sortedConsts[i][j]._coeff[0]));
		}
	}

	/* Copy velocity to the extended vector */
	for (int j = 0; j < N; ++j)
		_b[j + 1] = _bb[j];
}

bool ProfileHeuristic::test(double b0, double b1, double b2, int i)
{
	bool ok = false;
	if (abs(b1 - _velConsts[i]._upperLimit) < 10e-8)
		ok = true;
	else
	{
		for (size_t j = 0; j < _sortedConsts[i - 1].size(); ++j)
		{
			if (abs(b1 - (_sortedConsts[i - 1][j]._coeff[1] * b0 + _sortedConsts[i - 1][j]._upperLimit)) < 10e-8)
			{
				ok = true;
				break;
			}
		}

		for (size_t j = 0; j < _sortedConsts[i].size(); ++j)
		{
			if (abs(b1 - (b2 - _sortedConsts[i][j]._lowerLimit) / _sortedConsts[i][j]._coeff[1]) < 10e-8)
			{
				ok = true;
				break;
			}
		}

	}
	return ok;
}

void ProfileHeuristic::calculateTrajectoryImpl()
{
	//Acceleration calculation
	size_t N = _path.getLength();
	Path accel(_path.getDimension());
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
	_traj.SetAcceleration(accel);

	//Velocity calculation
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