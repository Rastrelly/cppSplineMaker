#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <cmath>
#include <ctime>
#include "TGaussSolver.h"
#include "TEigenSolver.h"
#include <GL/freeglut.h>

bool useEigen = true;


struct pt { double x = 0, y = 0; };
typedef std::vector<pt> ptVec;
typedef std::vector<std::vector<double>> splineCoeffVec;


ptVec oglBase = {}, oglSpline = {};
pt oglInterpPoint = {};

void renderChartVector(ptVec& vecRef, double r, double g, double b, int width)
{
	glColor3d(r, g, b);
	glLineWidth(width);

	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < vecRef.size(); ++i)
	{
		glVertex2d(
			vecRef[i].x,
			vecRef[i].y
		);
	}
	glEnd();
}


void drawCrossAtPoint(pt& crossPoint, double miX, double miY,
	double maX, double maY, double r, double g, double b,
	int width)
{
	glColor3d(r, g, b);
	glLineWidth(width);

	glBegin(GL_LINES);

	glVertex2d(miX, crossPoint.y);
	glVertex2d(maX, crossPoint.y);
	glVertex2d(crossPoint.x, miY);
	glVertex2d(crossPoint.x, maY);

	glEnd();
}

void cbDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	double xMin = 0;
	double xMax = 0;

	double yMin = 0;
	double yMax = 0;

	xMin = oglSpline.front().x;
	xMax = oglSpline.back().x;

	yMin = oglSpline[0].y;
	yMax = oglSpline[0].y;

	for (int i = 0; i < oglSpline.size(); ++i)
	{
		if (oglSpline[i].y < yMin) yMin = oglSpline[i].y;
		if (oglSpline[i].y > yMax) yMax = oglSpline[i].y;
	}

	double diap = yMax - yMin;

	yMin -= 0.5 * diap;
	yMax += 0.5 * diap;

	gluOrtho2D(xMin, xMax, yMin, yMax);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//draw charts
	renderChartVector(oglBase, 1, 0, 0, 1);
	renderChartVector(oglSpline, 0, 1, 1, 2);

	//draw interpolated point
	drawCrossAtPoint(oglInterpPoint, xMin, yMin, xMax, yMax, 0, 0, 1, 1);

	glutSwapBuffers();
}

void makeRandomVector(int vecSize, ptVec& outputVec)
{
	outputVec.clear();
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> xStepDist(1, 10);
	std::uniform_real_distribution<double> yDist(-100, 100);
	double x = yDist(gen);
	for (int i = 0; i < vecSize; ++i)
	{
		outputVec.push_back(
			{
				x,
				yDist(gen)
			}
		);
		x += xStepDist(gen);
	}
}

void makeCubSplineCoeffs(ptVec& inpDat, splineCoeffVec& coeffs)
{
	std::cout << "Making coeff search preps...\n";
	int l = inpDat.size();
	coeffs.clear();
	coeffs.resize(4);
	//a
	std::cout << "Starting calculating A...\n";
	for (int i = 0; i < l; ++i)
		coeffs[0].push_back(inpDat[i].y);
	//c
	std::cout << "Starting calculating C...\n";
	//need to make a system for c
	//get intermediate h
	std::cout << " - Defining H...\n";
	std::vector<double> h(l);
	h[0] = 0;
	for (int i = 1; i < l; ++i)
		h[i] = inpDat[i].x - inpDat[i - 1].x;
	//now make system for calculating c
	int sysSize = l - 2;
	splineCoeffVec cSetup(sysSize);

	std::cout << " - Building cSetup...\n";
	std::cout << " - cSetup size = " << sysSize << "\n";
	std::cout << " - H size = " << h.size() << "\n";
	for (int i = 0; i < sysSize; ++i)
	{
		cSetup[i].resize(sysSize + 1);
		for (int j = 0; j < cSetup[i].size(); ++j) cSetup[i][j] = 0;

		cSetup[i].back() = 3 * (((inpDat[i + 2].y - inpDat[i + 1].y) / h[i + 2]) -
			((inpDat[i + 1].y - inpDat[i].y) / h[i + 1]));

		std::cout << " -  - Section " << i << " at " << i + 1 << ", " << i + 2 << "\n";

		if (i - 1 >= 0)
			cSetup[i][i - 1] = h[i + 1];
		if (i < sysSize)
			cSetup[i][i] = 2 * (h[i + 1] + h[i + 2]);
		if (i + 1 < sysSize)
			cSetup[i][i + 1] = h[i + 2];
	}

	std::cout << " - cSetup built: \n";
	for (int i = 0; i < cSetup.size(); ++i)
	{
		std::cout << " - - " << i << "] ";
		for (int j = 0; j < cSetup[i].size(); ++j)
		{
			std::cout << "(" << j << ")" << cSetup[i][j] << " ";
		}
		std::cout << "\n";
	}

	std::cout << " - Running Gauss...\n";
	std::vector<double> cRoots = {};

	if (!useEigen)
	{
		TGaussSolver CFinder(cSetup);
		CFinder.solve();
		CFinder.returnRoots(cRoots);
	}
	else
	{
		TEigenSolver CFinder(cSetup);
		CFinder.returnRoots(cRoots);
	}

	std::cout << " - Roots found: ";
	for (double cr : cRoots)
		std::cout << cr << "; ";
	std::cout << "\n";

	coeffs[2].push_back(0);

	for (int i = 0; i < cRoots.size(); ++i)
	{
		coeffs[2].push_back(cRoots[i]);
	}

	coeffs[2].push_back(0);

	//d
	std::cout << "Starting calculating D...\n";
	coeffs[3].push_back(0);
	for (int i = 1; i < l; ++i)
	{
		coeffs[3].push_back((coeffs[2][i] - coeffs[2][i - 1]) / (3 * h[i]));
	}

	//b
	std::cout << "Starting calculating B...\n";
	coeffs[1].push_back(0);
	for (int i = 1; i < l; ++i)
	{
		coeffs[1].push_back
		(
			((inpDat[i].y - inpDat[i - 1].y) / h[i]) + coeffs[2][i] * h[i] - coeffs[3][i] * pow(h[i], 2)
		);
	}
	std::cout << "Coeff matrix complete!\n";
	for (int i = 0; i < coeffs.size(); ++i)
	{
		std::cout << i << "] ";
		for (int j = 0; j < coeffs[i].size(); ++j)
		{
			std::cout << "(" << j << ")" << coeffs[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void calcSplinePoints(ptVec& inpDat, ptVec& outpDat, int nPoints)
{
	splineCoeffVec coeffs = {};
	outpDat.clear();
	makeCubSplineCoeffs(inpDat, coeffs);
	double x1 = inpDat.front().x;
	double x2 = inpDat.back().x;
	double step = (x2 - x1) / (nPoints - 1);

	for (int i = 0; i < nPoints; ++i)
	{
		double cx = x1 + i * step;
		//find spline section index
		int ind = inpDat.size() - 1;
		for (int j = 1; j < inpDat.size(); ++j)
		{
			if (inpDat[j].x >= cx)
			{
				ind = j;
				break;
			}
		}
		double cy = coeffs[0][ind] +						//a
			coeffs[1][ind] * (cx - inpDat[ind].x) +			//b
			coeffs[2][ind] * pow(cx - inpDat[ind].x, 2) +	//c
			coeffs[3][ind] * pow(cx - inpDat[ind].x, 3);	//d	
		outpDat.push_back({ cx,cy });
	}
}

//calculate interpolated value
void getSplineInterpol(const ptVec& inpDat, const splineCoeffVec coeffs, double x, double &y)
{
	//find spline section index
	int ind = inpDat.size() - 1;
	for (int j = 1; j < inpDat.size(); ++j)
	{
		if (inpDat[j].x >= x)
		{
			ind = j;
			break;
		}
	}

	//get y
	y = coeffs[0][ind] +								//a
		coeffs[1][ind] * (x - inpDat[ind].x) +			//b
		coeffs[2][ind] * pow(x - inpDat[ind].x, 2) +	//c
		coeffs[3][ind] * pow(x - inpDat[ind].x, 3);		//d	
}

int main(int argc, char** argv)
{
	ptVec inputData = {};
	ptVec splineData = {};
	makeRandomVector(6, inputData);
	calcSplinePoints(inputData, splineData, 1000);
	std::cout << "+++ INPUT VECTOR +++\n";
	for (pt& cpt : inputData)
	{
		std::cout << cpt.x << "; " << cpt.y << ";\n";
	}

	std::cout << "+++ SPLINE VECTOR +++\n";
	for (pt& cpt : splineData)
	{
		std::cout << cpt.x << "; " << cpt.y << ";\n";
	}

	//get random interp value
	std::mt19937 gen(time(NULL));
	std::uniform_real_distribution<double> xRnd(inputData.front().x, 
												inputData.back().x);
	pt intPt = {0,0};
	intPt.x = xRnd(gen);
	
	splineCoeffVec coeffs = {};
	makeCubSplineCoeffs(inputData, coeffs);

	getSplineInterpol(inputData, coeffs, intPt.x, intPt.y);

	std::cout << "Interpolated: x = " << intPt.x << "; y = " << intPt.y << "\n";

	oglInterpPoint = intPt;

	std::cout << "+++ RENDERING +++\n";
	oglBase = inputData;
	oglSpline = splineData;

	glutInit(&argc, argv);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Chart");
	glClearColor(0.0, 0.0, 0.0, 0.0);

	glutDisplayFunc(cbDisplay);

	glutMainLoop();

}