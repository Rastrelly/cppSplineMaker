#pragma once

#include <iostream>
#include <vector>
#include <iterator>

class TGaussSolver
{
private:
	std::vector<std::vector<double>> A = {};
	std::vector<std::vector<double>> AVuln = {};
	std::vector<double> B = {};
	std::vector<double> BVuln = {};
	int n;
	std::vector<double> X = {};

	void prepSolver(std::vector<std::vector<double>> eqMatrix, bool showDebug);
	bool doDebugPrint = false;

	void debugPrintA();	
public:
	TGaussSolver(std::vector<std::vector<double>> eqMatrix);
	TGaussSolver(std::vector<std::vector<double>> eqMatrix, bool showDebug);
	void solve();
	void returnRoots(std::vector<double>& roots);
};

