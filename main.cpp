#include <Eigen/Core>
#include <limits>
#include <iostream>
#include <igl/signed_distance.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/copyleft/marching_cubes.h>

Eigen::VectorXd gridvals;
Eigen::MatrixXd P;
double isovalue;
int m;
igl::viewer::Viewer *viewer;

int valIndex(int i, int j, int k, int m)
{
    return i + j*m + k*m*m;
}

double cellGradientNormSquared(int x, int y, int z, double cubewidth, int m, const Eigen::VectorXd &gridvals)
{
    double ret = 0;
    double cellWidth = cubewidth / double(m);
    double area = cellWidth*cellWidth*cellWidth;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                ret += gridvals[valIndex(x + i, y + j, z + k, m)]*gridvals[valIndex(x + i, y + j, z + k, m)]* area / 3.0;
                int ni = (i + 1) % 2;
                int nj = (j + 1) % 2;
                int nk = (k + 1) % 2;
                ret += gridvals[valIndex(x + i, y + j, z + k, m)]*gridvals[valIndex(x + i, y + nj, z + nk, m)]* -area / 12.0;
                ret += gridvals[valIndex(x + i, y + j, z + k, m)]*gridvals[valIndex(x + ni, y + j, z + nk, m)]* -area / 12.0;
                ret += gridvals[valIndex(x + i, y + j, z + k, m)]*gridvals[valIndex(x + ni, y + nj, z + k, m)]* -area / 12.0;
                ret += gridvals[valIndex(x + i, y + j, z + k, m)]*gridvals[valIndex(x + ni, y + nj, z + nk, m)] * -area / 12.0;
            }
        }
    }
    return ret;
}

double sdfGradientResidual(double cubewidth, int m, const Eigen::VectorXd &gridvals)
{
    double ret = 0;
    for (int i = 0; i < m - 1; i++)
    {
        for (int j = 0; j < m - 1; j++)
        {
            for (int k = 0; k < m - 1; k++)
            {
                double grad = cellGradientNormSquared(i, j, k, cubewidth, m, gridvals);                
                double term = fabs(grad - 1.0);
                ret += term;
            }
        }
    }
    return ret;
}

void recomputeIsosurface()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXd vals = gridvals;
    for (int i = 0; i < vals.size(); i++)
        vals[i] -= isovalue;
    igl::copyleft::marching_cubes(vals, P, m, m, m, V, F);
    viewer->data.clear();
    viewer->data.set_mesh(V, F);
}

void fitIntoCube(Eigen::MatrixXd &V, double cubeWidth)
{
	double mins[3], maxs[3];
	for(int i=0; i<3; i++)
	{
		mins[i] = std::numeric_limits<double>::infinity();
		maxs[i] = -std::numeric_limits<double>::infinity();
	}
	int nverts = V.rows();
	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<3; j++)
		{
			mins[j] = std::min(mins[j], V(i,j));
			maxs[j] = std::max(maxs[j], V(i,j));
		}
	}
	double diameter = 0;
	for(int i=0; i<3; i++)
		diameter = std::max(diameter, maxs[i]-mins[i]);

	for(int i=0; i<nverts; i++)
	{
		for(int j=0; j<3; j++)
		{
			V(i,j) -= mins[j];
		}
	}
	
	double s = cubeWidth/diameter;
	if(diameter == 0.0)
		s = 1.0;
	V *= s;
}

void generateGridPoints(double cubeWidth, int m, Eigen::MatrixXd &P)
{
    P.resize(m*m*m, 3);
    for(int i=0; i<m; i++)
    {
        for(int j=0; j<m; j++)
        {
            for(int k=0; k<m; k++)
            {
                P(i + j*m + k*m*m, 0) = double(i)*cubeWidth/double(m-1);
                P(i + j*m + k*m*m, 1) = double(j)*cubeWidth/double(m-1);
                P(i + j*m + k*m*m, 2) = double(k)*cubeWidth/double(m-1);
            }
        }
    }
}

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cerr << "Usage: meshsdf [.txt file] [bounding cube width] [grid dimension m]" << std::endl;
		return -1;
	}
	

	double cubewidth = strtod(argv[2], NULL);
	m = strtol(argv[3], NULL, 10);

    std::ifstream ifs(argv[1]);
	if(!ifs)
	{
		std::cerr << "Couldn't read file " << argv[1] << std::endl;
		return -1;
	}

    gridvals.resize(m*m*m);
    for (int i = 0; i < m*m*m; i++)
    {
        ifs >> gridvals[i];
    }
	
    generateGridPoints(cubewidth, m, P);
	

    isovalue = 0;

    viewer = new igl::viewer::Viewer;

    recomputeIsosurface();

    double cubevol = double(m*m*m);

    std::cout << sdfGradientResidual(cubewidth, m, gridvals)/cubevol << std::endl;

    viewer->callback_init = [&](igl::viewer::Viewer& v)
    {
        v.ngui->addVariable("Isovalue", isovalue);
        v.ngui->addButton("Recompute Isosurface", recomputeIsosurface);
        v.screen->performLayout();
        return false;
    };

    viewer->data.set_face_based(true);
    viewer->launch();
}
