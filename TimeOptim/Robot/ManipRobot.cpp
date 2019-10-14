#include <algorithm>
#include "Robot\ManipRobot.h"
#include "Path\Path.h"

using namespace std;

Point ManipRobot::getAccelerationConstraint() const
{
	Point aMax;
	for (int i = 0; i < _dim; ++i)
		aMax.addCoord(_maxAccel);

	return aMax;
}

std::vector<double> ManipRobot::getAConstraintForPath(const Path& path) const
{
	vector<double> aMax;

	auto s_m_prime = path.ApproximateMiddleFirstDerivative();

	for (size_t i = 0; i < path.getLength() - 1; ++i)
	{
		double m = sqrt((s_m_prime[i] * s_m_prime[i]).getMaximum());
		double activeLimit = _maxAccel / m;
		aMax.push_back(activeLimit);
	}

	return aMax;
}

vector<double> ManipRobot::getVelocityConstraint() const
{
	vector<double> bMax;
	for (int i = 0; i < _dim; ++i)
		bMax.push_back(_maxVel);

	return bMax;
}

vector<double> ManipRobot::getBConstraintForPath(const Path& path) const
{
	vector<double> bMax;

	auto s_prime = path.ApproximateFirstDerivative();

	bMax.push_back(0.0); //Start velocity
	for (size_t i = 1; i < path.getLength() - 1; ++i)
	{
		double m = ((s_prime[i] * s_prime[i]).getMaximum());
		double activeLimit = _maxVel*_maxVel / m;
		bMax.push_back(activeLimit);
	}
	bMax.push_back(0.0); //End velocity

	return bMax;
}


Point ManipRobot::getTorgueConstraint() const
{
	Point bMax;
	for (int i = 0; i < _dim; ++i)
		bMax.addCoord(_maxTorgue);

	return bMax;
}

std::vector<Point> ManipRobot::getUConstraintForPath(const Path& path) const
{
	Point p;
	for (int i = 0; i < path.getDimension(); ++i)
		p.addCoord(_maxTorgue);
	return vector<Point>(path.getLength() - 2, p);
}

Point ManipRobot::GetDynM(const Point& s, const Point& sd) const
{
	auto M = GetInertiaMatrix(s);

	Point m;
	for (int j = 0; j < _dim; ++j)
	{
		double a = 0.0;
		for (int i = 0; i < _dim; ++i)
			a += M[j][i] * sd[i];
		m.addCoord(a);
	}
	return m;
}

Point ManipRobot::GetDynC(const Point& s, const Point& sd, const Point& sdd) const
{
	auto M = GetInertiaMatrix(s);
	auto C = GetCoriolisMatrix(s);

	Point c;
	for (int j = 0; j < _dim; ++j)
	{
		double a = 0.0;
		for (int i = 0; i < _dim; ++i)
		{
			a += M[j][i] * sdd[i];
			for (int k = 0; k < _dim; ++k)
			{
				a += C[j][i][k] * sd[k] * sd[i];
			}
		}
		c.addCoord(a);
	}
	return c;
}

Point ManipRobot::GetDynD(const Point& s) const
{
	return GetGravity(s);
}


std::vector<std::vector<double>> ManipRobot::GetInertiaMatrix(const Point& p) const
{
	double s2 = sin(p[1]);
	double c2 = cos(p[1]);
	double c3 = cos(p[2]);
	double s23 = sin(p[1] + p[2]);
	double c23 = cos(p[1] + p[2]);

	vector<vector<double>> mass;
	{
		vector<double> m;
		m.push_back(Iy2*s2*s2 + Iy3*s23*s23 + Iz1 + Iz2*c2*c2 + Iz3*c23*c23 + m2*r1*r1*c2*c2 + m3*(l1*c2 + r2*c23)*(l1*c2 + r2*c23));
		m.push_back(0.0);
		m.push_back(0.0);
		mass.push_back(m);
	}
	{
		vector<double> m;
		m.push_back(0.0);
		m.push_back(Ix2 + Ix3 + m3*l1*l1 + m2*r1*r1 + m3*r2*r2 + 2 * m3*l1*r2*c3);
		m.push_back(Ix3 + m3*r2*r2 + m3*l1*r2*c3);
		mass.push_back(m);
	}
	{
		vector<double> m;
		m.push_back(0.0);
		m.push_back(Ix3 + m3*r2*r2 + m3*l1*r2*c3);
		m.push_back(Ix3 + m3*r2*r2);
		mass.push_back(m);
	}
	return mass;
}

std::vector<std::vector<std::vector<double>>> ManipRobot::GetCoriolisMatrix(const Point& p) const
{
	vector<std::vector<std::vector<double>>> C;
	double s2 = sin(p[1]);
	double s3 = sin(p[2]);
	double c2 = cos(p[1]);
	double s23 = sin(p[1] + p[2]);
	double c23 = cos(p[1] + p[2]);
	//Init with zeros
	for (int i = 0; i < _dim; ++i)
	{
		vector<vector<double>> m;
		for (int j = 0; j < _dim; ++j)
		{
			vector<double> n;
			for (int k = 0; k < _dim; ++k)
				n.push_back(0.0);
			m.push_back(n);
		}
		C.push_back(m);
	}

	//111-113
	C[0][0][1] = (Iy2 - Iz2 - m2*r1*r1)*c2*s2 + (Iy3 - Iz3)*c23*s23 - m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
	C[0][0][2] = (Iy3 - Iz3)*c23*s23 - m3*r2*s23*(l1*c2 + r2*c23);
	//121-123 
	C[0][1][0] = (Iy2 - Iz2 - m2*r1*r1)*c2*s2 + (Iy3 - Iz3)*c23*s23 - m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
	//131-133
	C[0][2][0] = (Iy3 - Iz3)*c23*s23 - m3*r2*s23*(l1*c2 + r2*c23);
	//211-213
	C[1][0][0] = (Iz2 - Iy2 + m2*r1*r1)*c2*s2 + (Iz3 - Iy3)*c23*s23 + m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
	//221-223
	C[1][1][2] = -l1*m3*r2*s3;
	//231-233
	C[1][2][1] = -l1*m3*r2*s3;
	C[1][2][2] = -l1*m3*r2*s3;
	//311-313
	C[2][0][0] = (Iz3 - Iy3)*c23*s23 + m3*r2*s23*(l1*c2 + r2*c23);
	//321-323
	C[2][1][1] = l1*m3*r2*s3;

	return C;
}

Point ManipRobot::GetGravity(const Point& p) const
{
	Point g;

	double c2 = cos(p[1]);
	double c23 = cos(p[1] + p[2]);

	g.addCoord(0);
	g.addCoord(-(m2*9.81*r1 + m3*9.81*l1)*c2 - m3*9.81*r2*c23);
	g.addCoord(-m3*9.81*r2*c23);
		
	return g;
}

std::vector<double> ManipRobot::getParams() const
{
	vector<double> a;
	//Ix
	a.push_back(Ix1);
	a.push_back(Ix2);
	a.push_back(Ix3);
	//Iy
	a.push_back(Iy1);
	a.push_back(Iy2);
	a.push_back(Iy3);
	//Iz
	a.push_back(Iz1);
	a.push_back(Iz2);
	a.push_back(Iz3);
	//m
	a.push_back(m1);
	a.push_back(m2);
	a.push_back(m3);
	//r
	a.push_back(r0);
	a.push_back(r1);
	a.push_back(r2);
	//l
	a.push_back(l0);
	a.push_back(l1);
	a.push_back(l2);

	return a;
}