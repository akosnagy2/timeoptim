#include "Profile\ProfileGurobi.h"
#include "Profile\ProfileHeuristic.h"
#include "Robot\ManipRobot.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void cycle(Profile& prof, std::ostream &log, bool useDynamics)
{
	log << "------------------" << endl;
	log << "Profile: " << prof.getName() << endl;

	if (prof.createModel(useDynamics) && prof.solve())
	{
		log << prof.getName() << " Model creation sucessfull, runtime: " << prof.getModelCreationTimes().back() * 1000.0 << "ms" << endl;
		log << prof.getName() << " Solver sucessfull, runtime: " << prof.getSolveTimes().back() * 1000.0 << "ms" << endl;		
	}
	else
	{
		log << "Solver failed." << endl;
	}
}

void check(const Profile& prof1, const Profile& prof2)
{
	const Trajectory& t1 = prof1.getTrajectory();
	const Trajectory& t2 = prof2.getTrajectory();

	if (t1.GetVelocity().getLength() != t2.GetVelocity().getLength())
	{
		cout << "Velocity-Path length different" << endl;
		return;
	}
	
	for (size_t i = 0; i < t1.GetPath().getLength(); ++i)
	{
		for (size_t j = 0; j < t1.GetPath().getDimension(); ++j)
		{
			if (!isfinite(t1.GetVelocity()[i][j]))
				cout << prof1.getName() << ": infinite solution at " << i << ", " << j << endl;
			if (!isfinite(t2.GetVelocity()[i][j]))
				cout << prof2.getName() << ": infinite solution at " << i << ", " << j << endl;
		}
		if (abs(prof1.getB()[i] - prof2.getB()[i]) > 10e-8)
			cout << prof1.getName() + "-" + prof2.getName() + " difference found at " << i << ", diff: " << abs(prof1.getB()[i] - prof2.getB()[i]) << endl;
	}
}

void main()
{
	try 
	{
		int N = 2000;
		string dir = "../Measurement/3DoF_Path/";
		string logDir = dir + "log/";
		string pathName = "path_" + std::to_string(N) + ".txt";
		Path path(dir + pathName);

		ManipRobot rob(100.5, 10.5, 100.5, path.getDimension());
		rob.l0 = 2.0000e-01;
		rob.l1 = 4.9126e-01 - rob.l0;
		rob.l2 = 1.9900e-01 + 1.2500e-01;
		rob.m1 = 1.5;
		rob.m2 = 1.2;
		rob.m3 = 1.0;
		rob.r0 = rob.l0 * 0.4;
		rob.r1 = rob.l1 * 0.4;
		rob.r2 = rob.l2 * 0.4;
		rob.Ix1 = 5 * rob.m1;
		rob.Iy1 = 5 * rob.m1;
		rob.Iz1 = 5 * rob.m1;
		rob.Ix2 = 4.75 * rob.m2;
		rob.Iy2 = 4.75 * rob.m2;
		rob.Iz2 = 4.75 * rob.m2;
		rob.Ix3 = 4.75 * rob.m3;
		rob.Iy3 = 4.75 * rob.m3;
		rob.Iz3 = 4.75 * rob.m3;

		ProfileHeuristic prof_SLP(path, rob);
		cycle(prof_SLP, std::cout, true);

		ProfileGurobi prof_SOCP(path, rob, false);
		cycle(prof_SOCP, std::cout, true);

		ProfileGurobi prof_LP(path, rob, true);
		cycle(prof_LP, std::cout, true);

		check(prof_LP, prof_SOCP);
		check(prof_SLP, prof_LP);

		const Trajectory& traj = prof_SLP.getTrajectory();
		traj.GetPath().savePath(logDir + "path_sim", traj.GetTime());
		traj.GetVelocity().savePath(logDir + "vel_sim", traj.GetTime());
		traj.GetAcceleration().savePath(logDir + "vel_sim", traj.GetTime());
		traj.GetDynamics().savePath(logDir + "vel_sim", traj.GetTime());
	}
	catch (exception e)
	{
		cout << "Program failed: " << e.what() << endl;
	}
}