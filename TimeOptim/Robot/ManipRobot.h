#pragma once
#include "Robot.h"

class ManipRobot : public Robot
{
public:
	ManipRobot(double maxVel, double maxAccel, double maxTorgue, int dim)
		: _maxVel(maxVel), _maxAccel(maxAccel), _maxTorgue(maxTorgue), Robot(dim) {};
	Point getAccelerationConstraint() const;
	std::vector<double> getAConstraintForPath(const Path& path) const;
	std::vector<double> getVelocityConstraint() const;
	std::vector<double> getBConstraintForPath(const Path& path) const;
	Point getTorgueConstraint() const;
	std::vector<Point> getUConstraintForPath(const Path& path) const;
	Point GetDynM(const Point& s, const Point& sd) const;
	Point GetDynC(const Point& s, const Point& sd, const Point& sdd) const;
	Point GetDynD(const Point& s) const;
	std::vector<std::vector<double>> GetInertiaMatrix(const Point& p) const;
	std::vector<std::vector<std::vector<double>>> GetCoriolisMatrix(const Point& p) const;
	Point GetGravity(const Point& p) const;
	std::vector<double> getParams() const;
public:
	double _maxVel, _maxAccel, _maxTorgue;
	double Ix1, Ix2, Ix3, Iy1, Iy2, Iy3, Iz1, Iz2, Iz3;
	double m1, m2, m3;
	double r0, r1, r2;
	double l0, l1, l2;
};