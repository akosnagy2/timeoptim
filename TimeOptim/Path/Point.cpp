#include "Point.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <utility>
#include <algorithm>
#include <exception>

Point::Point(std::vector<double> p) {
	_point = p;
}

Point::Point(int dim, double val) : Point(std::vector<double>(dim, val))
{
}

void Point::addCoord(double p) {
	_point.push_back(p);
}

void Point::removeCoord(int index)
{
	_point.erase(_point.begin() + index);
}

void Point::removeLastCoord()
{
	removeCoord(_point.size() - 1);
}

bool Point::IsDimensionEqual(const Point& b) const
{
	if (b.getDimension() != this->getDimension())
		return false;
	else
		return true;
}

const Point Point::operator+(const Point& b) const
{
	Point r;

	if (!IsDimensionEqual(b))
		throw std::exception("Dimension mismatch.");

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) + b.getCoord(i));

	return r;
}

const Point Point::operator-(const Point& b) const
{
	Point r;

	if (!IsDimensionEqual(b))
		throw std::exception("Dimension mismatch.");

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) - b.getCoord(i));

	return r;
}

const Point Point::operator*(const Point& b) const
{
	Point r;

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) * b[i]);

	return r;
}

const Point Point::operator/(const Point& b) const
{
	Point r;

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) / b[i]);

	return r;
}

const Point Point::operator*(double b) const
{
	Point r;

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) * b);

	return r;
}

const Point Point::operator/(double b) const
{
	Point r;

	for (int i = 0; i < this->getDimension(); ++i)
		r.addCoord(this->getCoord(i) / b);

	return r;
}

const Point Point::Sqrt(const Point& b)
{
	Point a;
	for (int j = 0; j < b.getDimension(); ++j)
		a.addCoord(sqrt(b[j]));

	return a;
}

double Point::getMaximum() const
{
	return *std::max_element(_point.begin(), _point.end());
}

double Point::getMinimum() const
{
	return *std::min_element(_point.begin(), _point.end());
}