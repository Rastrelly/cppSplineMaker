#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class TEigenSolver
{
private:	
	Eigen::VectorXd x;
public:

	//constructs and solves the system
	TEigenSolver(std::vector<std::vector<double>> eqData)
	{
		//init containers
		Eigen::MatrixXd A(eqData.size(), eqData.front().size()-1);
		Eigen::VectorXd B(eqData.size());
		//fill vector B
		int cnt = 0;
		for (std::vector<double>& v : eqData)
		{			
			B[cnt] = v.back();
			cnt++;
		}
		//fill matrix A
		for (int i = 0; i < eqData.size(); ++i)
		{
			for (int j = 0; j < eqData.front().size() - 1; ++j)
			{
				A(i, j) = eqData[i][j];
			}
		}
		//run solver
		auto lu = A.fullPivLu();
		x = lu.solve(B);
	}

	//returns roots vector
	void returnRoots(std::vector<double>& roots)
	{
		roots.clear();
		int c = x.size();
		for (int j = 0; j < c; ++j)
		{
			if (abs(x(j)) < 0.00000001)
			{
				roots.push_back(0);
				continue;
			}
			roots.push_back(x(j));
		}
	}

};

