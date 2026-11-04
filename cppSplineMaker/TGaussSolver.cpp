#include "TGaussSolver.h"


TGaussSolver::TGaussSolver(std::vector<std::vector<double>> eqMatrix)
{
	prepSolver(eqMatrix, false);
}

TGaussSolver::TGaussSolver(std::vector<std::vector<double>> eqMatrix, bool showDebug)
{
	prepSolver(eqMatrix, showDebug);
}

void TGaussSolver::prepSolver(std::vector<std::vector<double>> eqMatrix,
	bool showDebug)
{
	//assume matrix is okay
	doDebugPrint = showDebug;
	n = eqMatrix.size();
	A.clear();
	AVuln.clear();
	B.clear();
	BVuln.clear();
	X.clear();
	for (int i = 0; i < n; ++i)
	{
		A.push_back({});
		for (int j = 0; j < n; ++j)
		{
			A[i].push_back(eqMatrix[i][j]);
		}
		B.push_back(eqMatrix[i].back());
	}
	//we now have A, B, n and can run solver
	//copy A to vulnerable version
	std::copy(A.begin(), A.end(), std::back_inserter(AVuln));
	std::copy(B.begin(), B.end(), std::back_inserter(BVuln));
	X.resize(n);

	if (doDebugPrint)
	{
		std::cout << "Initial A matrix:\n";
		debugPrintA();
	}
}

void TGaussSolver::solve()
{
	for (int k = 0; k < n; ++k)
	{
		for (int j = k + 1; j < n; ++j)
		{
			AVuln[k][j] = AVuln[k][j] / AVuln[k][k];
		}
		BVuln[k] = BVuln[k] / AVuln[k][k];
		for (int i = k + 1; i < n; ++i)
		{
			for (int j = k + 1; j < n; ++j)
			{
				AVuln[i][j] -= AVuln[i][k] * AVuln[k][j];
			}
			BVuln[i] -= AVuln[i][k] * BVuln[k];
		}
		if (doDebugPrint)
			std::cout << "Iteration k = " << k << "\n";
		debugPrintA();
	}
	X[n - 1] = BVuln[n - 1] / AVuln[n - 1][n - 1];
	for (int i = n - 1; i >= 0; --i)
	{
		double S = 0;
		for (int j = i + 1; j < n; ++j)
		{
			S += AVuln[i][j] * X[j];
		}
		X[i] = BVuln[i] - S;
	}
}

void TGaussSolver::returnRoots(std::vector<double>& roots)
{
	std::copy(X.begin(), X.end(), std::back_inserter(roots));
}

void TGaussSolver::debugPrintA()
{
	if (doDebugPrint)
	{
		int i = 0;
		for (std::vector<double> row : AVuln)
		{
			for (double rval : row)
			{
				std::cout << rval << " ";
			}
			std::cout << BVuln[i];
			std::cout << "\n";
			i++;
		}
	}
}