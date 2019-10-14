#include "Trajectory.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

Trajectory::Trajectory(std::string fileName) :
	_path(0),
	_vel(0),
	_accel(0),
	_tau(0)
{
	string line;
	ifstream fs(fileName);

	cout << "------------------" << endl;
	cout << "Load trajectory from: " << fileName << endl;

	if (!fs.good())
	{
		cout << "Trajectory load failed: File not found." << endl;
		throw exception("Trajectory File not found");
	}
		
	bool first = true;
	while (getline(fs, line))
	{
		std::stringstream ss(line);
		Point pos;
		string a;
		while (ss.good())
		{
			ss >> a;
			pos.addCoord(stof(a));
		}
		if (first) //First point
		{
			_path = Path(pos.getDimension() - 1);
			_vel = Path(_path.getDimension());
			_accel = Path(_path.getDimension());
			_tau = Path(_path.getDimension());
			first = false;
		}

		double time = pos.getCoord(pos.getDimension() - 1);
		pos.removeLastCoord();
		AddPoint(pos, time);
	}

	if (_path.getLength() == 0)
	{
		cout << "Trajectory load failed: No path point found." << endl;
		throw exception("No path point found in file");
	}

	cout << "Trajectory loaded: " << _path.getLength() << " point, " << _path.getDimension() << " dimension" << endl;
}

void Trajectory::Clear()
{
	_path.clear();
	_time.clear();
	_vel.clear();
	_accel.clear();
	_tau.clear();
}

void Trajectory::AddPoint(const Point& point, double time)
{
	AddPoint(point, Point(point.getDimension(), 0.0), time);
}
void Trajectory::AddPoint(const Point& point, const Point& vel, double time)
{
	AddPoint(point, vel, Point(point.getDimension(), 0.0), time);
}
void Trajectory::AddPoint(const Point& point, const Point& vel, const Point& accel, double time)
{
	AddPoint(point, vel, accel, Point(point.getDimension(), 0.0), time);
}
void Trajectory::AddPoint(const Point& point, const Point& vel, const Point& accel, const Point& tau, double time)
{
	if (_path.getDimension() != point.getDimension())
		throw std::exception("Path point dimension mismatch.");

	if (_path.getDimension() != vel.getDimension())
		throw std::exception("Velocity point dimension mismatch.");

	if (_path.getDimension() != accel.getDimension())
		throw std::exception("Acceleration point dimension mismatch.");

	if (_path.getDimension() != tau.getDimension())
		throw std::exception("Dynamics point dimension mismatch.");

	_time.push_back(time);
	_path.AddPoint(point);
	_vel.AddPoint(vel);
	_accel.AddPoint(vel);
	_tau.AddPoint(vel);
}

void Trajectory::SetPath(const Path& path)
{
	if (_path.getDimension() != path.getDimension())
		throw std::exception("Path point dimension mismatch.");

	_path = path;

}
void Trajectory::SetVelocity(const Path& vel)
{
	if (_vel.getDimension() != vel.getDimension())
		throw std::exception("Velocity point dimension mismatch.");

	_vel = vel;
}
void Trajectory::SetAcceleration(const Path& accel)
{
	if (_accel.getDimension() != accel.getDimension())
		throw std::exception("Acceleration point dimension mismatch.");

	_accel = accel;
}
void Trajectory::SetDynamics(const Path& tau)
{
	if (_tau.getDimension() != tau.getDimension())
		throw std::exception("Dynamics point dimension mismatch.");

	_tau = tau;
}
void Trajectory::SetTime(const std::vector<double>& time)
{
	_time = time;
}

std::vector<double> Trajectory::GetCenterTime() const
{
	size_t n =  _time.size();
	vector<double> t;
	t.reserve(n - 1);

	for (int i = 0; i < n - 1; ++i)
		t.push_back((_time[i] + _time[i + 1]) / 2.0);

	return t;
}