#include "Path.h"
#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>  

using namespace std;

Path::Path(std::string fileName)
{
	string line;
	ifstream fs(fileName);

	cout << "------------------" << endl;
	cout << "Load path from: " << fileName << endl;

	if (!fs.good())
	{
		cout << "Path load failed: File not found." << endl;
		throw exception("Path File not found");
	}

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
		if (_path.size() == 0) //First point
			_dim = pos.getDimension();

		AddPoint(pos);
	}

	if (_path.size() == 0)
	{
		cout << "Path load failed: No path point found." << endl;
		throw exception("No path point found in file");
	}

	setEquidistantTheta();

	cout << "Path loaded: " << getLength() << " point, " << getDimension() << " dimension" << endl;
}

void Path::AddPoint(Point p)
{
	if (p.getDimension() != _dim)
		throw std::exception("Dimension mismatch.");

	_path.push_back(p);
}

Path Path::ApproximateMiddleFirstDerivative() const
{
	Path pp(_dim);
	for (size_t i = 1; i < getLength(); ++i)
	{
		Point e = (_path[i] - _path[i - 1]);		
		pp.AddPoint(e / getTheta());
	}

	return pp;
}

Path Path::ApproximateMiddleSecondDerivative() const
{
	Path pp(_dim);

	for (size_t i = 0; i < getLength() - 1; ++i)
	{
		Point e;

		if (i == 0)
			e = (_path[i] - _path[i + 1] * 2 + _path[i + 2]);
		else if (i == (getLength() - 2))
			e = (_path[i - 1] - _path[i] * 2 + _path[i + 1]);
		else if ((i == 1) || (i == (getLength() - 3)))
			e = (_path[i - 1] - _path[i] - _path[i + 1] + _path[i + 2]);
		else
			e = (_path[i - 2] * (-10.0) / 48.0 + _path[i - 1] * 26.0 / 16.0 - _path[i] * 34.0 / 24.0 - _path[i + 1] * 34.0 / 24.0 + _path[i + 2] * 26.0 / 16.0 - _path[i + 3] * 10.0 / 48.0);		

		pp.AddPoint(e / (2 * getTheta()*getTheta()));
	}

	return pp;
}

Path Path::ApproximateFirstDerivative() const
{
	Path pp(_dim);

	if (getLength() == 3)
	{
		pp.AddPoint((_path[0] * (-3.0) / 2.0 + _path[1] * 2.0 - _path[2] * 1.0 / 2.0) / getTheta());
		pp.AddPoint((_path[0] * (-1.0) / 2.0 + _path[2] * 1.0 / 2.0) / getTheta());
		pp.AddPoint((_path[0] * 1.0 / 2.0 - _path[1] * 2.0 + _path[2] * 3.0 / 2.0) / getTheta());
	}
	else if (getLength() == 4)
	{
		pp.AddPoint((_path[0] * (-3.0) / 2.0 + _path[1] * 2.0 - _path[2] * 1.0 / 2.0) / getTheta());
		pp.AddPoint((_path[0] * (-1.0) / 2.0 + _path[2] * 1.0 / 2.0) / getTheta());
		pp.AddPoint((_path[1] * (-1.0) / 2.0 + _path[3] * 1.0 / 2.0) / getTheta());
		pp.AddPoint((_path[0] * 1.0 / 2.0 - _path[1] * 2.0 + _path[2] * 3.0 / 2.0) / getTheta());
	}
	else
	{
		for (int i = 0; i < 2; ++i)
		{
			Point e = (_path[i] * (-11.0) / 6.0 + _path[i + 1] * 3.0 - _path[i + 2] * 3.0 / 2.0 + _path[i + 3] * 1.0 / 3.0);
			pp.AddPoint(e / getTheta());
		}

		for (size_t i = 2; i < getLength() - 2; ++i)
		{
			Point e = (_path[i - 2] * 1.0 / 12.0 - _path[i - 1] * 2.0 / 3.0 + _path[i + 1] * 2.0 / 3.0 - _path[i + 2] * 1.0 / 12.0);
			pp.AddPoint(e / getTheta());
		}

		for (size_t i = getLength() - 2; i < getLength(); ++i)
		{
			Point e = (_path[i - 2] * 1.0 / 2.0 - _path[i - 1] * 2.0 + _path[i] * 3.0 / 2.0);
			pp.AddPoint(e / getTheta());
		}
	}
	return pp;
}

Path Path::ApproximateSecondDerivative() const
{
	Path pp(_dim);

	if (getLength() == 3)
	{
		pp.AddPoint((_path[0] - _path[1] * 2.0 + _path[2]) / (getTheta()*getTheta()));
		pp.AddPoint((_path[0] - _path[1] * 2.0 + _path[2]) / (getTheta()*getTheta()));
		pp.AddPoint((_path[0] - _path[1] * 2.0 + _path[2]) / (getTheta()*getTheta()));
	}
	else
	{
		for (int i = 0; i < 2; ++i)
		{
			Point e = (_path[i] * (35.0) / 12.0 + _path[i + 1] * -(26.0 / 3.0) + _path[i + 2] * 19.0 / 2.0 - _path[i + 3] * 14.0 / 3.0 + _path[i + 4] * 11.0 / 12.0);
			pp.AddPoint(e / (getTheta()*getTheta()));
		}

		for (size_t i = 2; i < getLength() - 2; ++i)
		{
			Point e = (_path[i - 2] * (-1.0 / 12.0) + _path[i - 1] * 4.0 / 3.0 - _path[i] * 5.0 / 2.0 + _path[i + 1] * 4.0 / 3.0 - _path[i + 2] * 1.0 / 12.0);
			pp.AddPoint(e / (getTheta()*getTheta()));
		}

		for (size_t i = getLength() - 2; i < getLength(); ++i)
		{
			Point e = (_path[i - 2] * 1.0 - _path[i - 1] * 2.0 + _path[i] * 1.0);
			pp.AddPoint(e / (getTheta()*getTheta()));
		}
	}
	return pp;
}

Path Path::GetMiddlePath() const
{
	Path pp(_dim);
	for (size_t i = 1; i < getLength(); ++i)
	{
		Point e = (_path[i] + _path[i - 1]);
		pp.AddPoint(e / 2.0);
	}

	return pp;
}

void Path::savePath(std::string filePrefix) const
{
	std::ofstream file;
	file.open(filePrefix + ".txt");

	for (size_t i = 0; i < getLength(); i++)
	{
		for (int j = 0; j < _path[0].getDimension(); ++j)
			file << _path[i][j] << " ";
		file << std::endl;
	}
	file.close();
}

void Path::savePath(std::string filePrefix, std::vector<double> time) const
{
	std::ofstream file;
	file.open(filePrefix + ".txt");
	
	if (time.size() != getLength())
		return;

	for (size_t i = 0; i < getLength(); i++)
	{
		for (int j = 0; j < _path[0].getDimension(); ++j)
			file << _path[i][j] << " ";
		file << time[i];
		file << std::endl;
	}
	file.close();
}

std::vector<double> Path::getVector() const {
	std::vector<double> pp;
	for (size_t i = 0; i < getLength(); ++i)
		for (int j = 0; j < getDimension(); ++j)
			pp.push_back(_path[i][j]);
	return pp;
}

void Path::removePoint(const int i)
{
	assert(i >= 0);
	assert(i < _path.size());

	_path.erase(_path.begin() + i);
}