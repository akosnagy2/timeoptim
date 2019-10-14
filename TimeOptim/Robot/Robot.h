#pragma once
#include "Path\Path.h"

class Robot
{
public:
	Robot(int dimension) : _dim(dimension) {}
	int getDimension() const {
		return _dim;
	}
	virtual Point getAccelerationConstraint() const = 0;
	virtual std::vector<double> getAConstraintForPath(const Path& path) const = 0;
	virtual std::vector<double> getVelocityConstraint() const = 0;
	virtual std::vector<double> getBConstraintForPath(const Path& path) const = 0;
	virtual Point getTorgueConstraint() const = 0;
	virtual std::vector<double> getParams() const = 0;
	virtual std::vector<Point> getUConstraintForPath(const Path& path) const = 0;
	virtual Point GetDynM(const Point& s, const Point& sd) const = 0;
	virtual Point GetDynC(const Point& s, const Point& sd, const Point& sdd) const = 0;
	virtual Point GetDynD(const Point& s) const = 0;
protected:
	int _dim;
};