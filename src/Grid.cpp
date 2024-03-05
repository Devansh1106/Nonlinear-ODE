#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include "../include/Master.hpp"
using namespace std;

Grid::Grid(double xMinimum, double xMaximum, int mNodes)
{
	xMin = xMinimum; xMax = xMaximum; Nodes = mNodes;
	std::vector<double>::size_type zNumNodes = Nodes;
	double mpNodes = (double)(Nodes);
	double stepsize = (xMax-xMin)/(mpNodes-1.0);
	for (int i = 0; i < Nodes; i++)
	{
		Mesh.push_back(xMin + (i * stepsize));
	}assert(Mesh.size() == zNumNodes);
	// cout<<Mesh[Nodes-1];
}