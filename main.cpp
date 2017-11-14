#include <Eigen/Core>
#include <limits>
#include <iostream>
#include <igl/signed_distance.h>
#include <igl/readOBJ.h>

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

void computeSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double cubeWidth, int m, Eigen::VectorXd &gridvals)
{
	Eigen::MatrixXd P(m*m*m, 3);
	for(int i=0; i<m; i++)
	{
		for(int j=0; j<m; j++)
		{
			for(int k=0; k<m; k++)
			{
				P(i*m*m + j*m + k, 0) = double(i)*cubeWidth/double(m-1);
				P(i*m*m + j*m + k, 1) = double(j)*cubeWidth/double(m-1);
				P(i*m*m + j*m + k, 2) = double(k)*cubeWidth/double(m-1);
			}
		}
	}
	Eigen::VectorXi I;
	Eigen::MatrixXd C;
	Eigen::MatrixXd N;
	igl::signed_distance(P, V, F, igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), gridvals, I, C, N);
}

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cerr << "Usage: meshsdf [.obj file] [bounding cube width] [grid dimension m]" << std::endl;
		return -1;
	}
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::readOBJ(argv[1], V, F);

	double cubewidth = strtod(argv[2], NULL);
	int m = strtol(argv[3], NULL, 10);

	if(V.rows() == 0)
	{
		std::cerr << "Couldn't read OBJ file " << argv[1] << std::endl;
		return -1;
	}
	fitIntoCube(V, cubewidth);
	Eigen::VectorXd gridvals;
	computeSDF(V, F, cubewidth, m, gridvals);
	std::cout << gridvals << std::endl;

}
