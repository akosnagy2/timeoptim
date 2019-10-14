#pragma once
#include "Path.h"

class Trajectory
{
public:
	Trajectory(int dim) : 
		_path(dim), 
		_vel(dim), 
		_accel(dim), 
		_tau(dim) {}
	Trajectory(Path& p, std::vector<double> t) : 
		_path(p), 
		_vel(p.getDimension()), 
		_accel(p.getDimension()), 
		_tau(p.getDimension()) {}
	Trajectory(std::string fileName);

	int GetDimension() const { return _path.getDimension(); }
	void Clear();
	void AddPoint(const Point& point, double time);
	void AddPoint(const Point& point, const Point& vel, double time);
	void AddPoint(const Point& point, const Point& vel, const Point& accel, double time);
	void AddPoint(const Point& point, const Point& vel, const Point& accel, const Point& tau, double time);
	void SetPath(const Path& path);
	void SetVelocity(const Path& vel);
	void SetAcceleration(const Path& accel);
	void SetDynamics(const Path& tau);
	void SetTime(const std::vector<double>& time);
	const Path& GetPath() const { return _path; }
	const Path& GetVelocity() const { return _vel; }
	const Path& GetAcceleration() const { return _accel; }
	const Path& GetDynamics() const { return _tau; }
	const std::vector<double>& GetTime() const { return _time; }
	std::vector<double> GetCenterTime() const;
private:
	Path _path, _vel, _accel, _tau;
	std::vector<double> _time;
};