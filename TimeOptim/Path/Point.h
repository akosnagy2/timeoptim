#pragma once
#include <vector>

class Point
{
public:
	Point()  {}
	Point(std::vector<double> p);
	Point(int dim, double val);
	void addCoord(double p);
	void removeCoord(int index);
	void removeLastCoord();
	int getDimension() const {
		return (int)_point.size();
	}
	double getCoord(int index) const {
		return _point[index];
	}
	double getMaximum() const;
	double getMinimum() const;
	operator std::vector<double>() const {
		return _point;
	}
	virtual const Point operator+(const Point &b) const;
	virtual const Point operator-(const Point &b) const;
	virtual const Point	operator*(const Point& b) const;
	virtual const Point	operator/(const Point& b) const;
	virtual const Point operator*(double b) const;
	virtual const Point operator/(double b) const;
	double& operator[] (const int i) {
		return _point[i];
	}
	double operator[] (const int i) const {
		return _point[i];
	}
	static const Point Sqrt(const Point& b);	
	bool IsDimensionEqual(const Point& b) const;
protected:	
	std::vector<double> _point;
};