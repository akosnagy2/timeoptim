#pragma once

#include "Point.h"
#include <vector>

class Path
{
public:
	Path(int dim) : _dim(dim), _theta(0.0) {}
	Path(std::vector<Point> p) : _path(p), _theta(1.0 / (p.size() - 1)), _dim(p[0].getDimension()) {}
	Path(std::string fileName);
	virtual void AddPoint(Point p);
	void clear() {
		_path.clear();
	}
	Path ApproximateMiddleFirstDerivative() const;
	Path ApproximateMiddleSecondDerivative() const;
	Path ApproximateFirstDerivative() const;
	Path ApproximateSecondDerivative() const;
	Path GetMiddlePath() const;
	double getTheta() const {
		return _theta;
	}
	void setTheta(double t) { 
		_theta = t; 
	}
	void setEquidistantTheta() {
		_theta = 1.0 / (getLength() - 1.0);
	}
	size_t getLength() const { 
		return _path.size(); 
	}
	const Point& back() const {
		return _path.back();
	}
	int getDimension() const {
		return _dim;
	}
	std::vector<double> getVector() const;
	void removePoint(const int i);
	operator std::vector<Point>& () {
		return _path;
	}
	Point& operator[] (const int i) {
		return _path[i];
	}
	const Point& operator[] (const int i) const {
		return _path[i];
	}
	void savePath(std::string filePrefix) const;
	void savePath(std::string filePrefix, std::vector<double> time) const;
protected:
	int _dim;
	double _theta;
	std::vector<Point> _path;
};