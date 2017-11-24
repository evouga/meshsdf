#include <Eigen/Core>
#include <limits>
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>
#include <igl/copyleft/marching_cubes.h>
#include "SDFAlgs.h"

Eigen::VectorXd gridvals;
Eigen::MatrixXd P;
double isovalue;
int m;
igl::viewer::Viewer *viewer;

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


struct options
{
    bool skipgen;
    bool outputsdf;
    bool showviewer;
    double cubeWidth;
    double inflateFactor;
    double n[3];
    double u[3];
};

bool parseFlags(int argc, char *argv[], options &opt)
{
    opt.skipgen = false;
    opt.outputsdf = false;
    opt.showviewer = false;
    opt.cubeWidth = 1.0;
    opt.inflateFactor = std::numeric_limits<double>::infinity();
    opt.n[0] = 0;
    opt.n[1] = 0;
    opt.n[2] = -1;
    opt.u[0] = 0;
    opt.u[1] = 1;
    opt.u[2] = 0;
    for (int arg = 3; arg < argc; arg++)
    {
        if (argv[arg][0] != '-')
            return false;
        switch (argv[arg][1])
        {
        case 's':
            opt.skipgen = true;
            break;
        case 'o':
            opt.outputsdf = true;
            break;
        case 'v':
            opt.showviewer = true;
            break;
        case 'w':
            if (arg + 1 >= argc)
                return false;
            arg++;
            opt.cubeWidth = strtod(argv[arg], NULL);
            break;
        case 't':
            if (arg + 1 >= argc)
                return false;
            arg++;
            opt.inflateFactor = strtod(argv[arg], NULL);
            break;
        case 'c':
            if (arg + 6 >= argc)
                return false;
            arg++;
            for (int i = 0; i < 3; i++)
            {
                opt.n[i] = strtod(argv[arg], NULL);
                arg++;
            }
            for (int i = 0; i < 3; i++)
            {
                opt.u[i] = strtod(argv[arg], NULL);
                arg++;
            }
            break;
        default:
            return false;
        }
    }
    return true;
}

void printUsage()
{
    std::cerr << "Usage: meshsdf (input-file) (grid dimension m) [flags]" << std::endl;
    std::cerr << "Computes a signed distance field from the mesh file input-file and optionally visualizes its zero level set." << std::endl;
    std::cerr << "Scales the model to fit into a unit cube, then builds an m x m x m grid for computing the SDF." << std::endl;
    std::cerr << "Flags can be:" << std::endl;
    std::cerr << " -s:  Skip SDF generation. In this case the input file should be a text file SDF instead of a mesh." << std::endl;
    std::cerr << " -o:  Print the SDF to the consoele." << std::endl;
    std::cerr << " -v:  Launch the visualization of the SDF after generation." << std::endl;
    std::cerr << " -w (width):  Scales the object to fit a bounding cube of the given side length (default 1.0)." << std::endl;
    std::cerr << " -t (factor):  Thickens the object by the given fraction of the bounding cube length before generating the SDF." << std::endl;
    std::cerr << " -c (nx) (ny) (nz) (ux) (uy) (uz): Generate a signed distance field as if viewing the object from a camera whose \"in\" direction is (nx, ny, nz) (default: (0,0,-1)) and whose \"up\" direction is (ux, uy, uz) (default: (0,1,0))." << std::endl;
}

Eigen::Matrix3d cameraMatrix(const options &opt)
{
    Eigen::Vector3d n(opt.n[0], opt.n[1], opt.n[2]);
    Eigen::Vector3d u(opt.u[0], opt.u[1], opt.u[2]);
    Eigen::Vector3d r = n.cross(u);
    double rnorm = r.norm();
    if (rnorm < 1e-6)
    {
        std::cerr << "Error: camera n and u directions degenerate. Using default settings." << std::endl;
        Eigen::Matrix3d ret;
        ret.setIdentity();
        return ret;
    }

    r /= rnorm;
    n /= n.norm();    

    if (n.dot(u) > 1e-6)
    {
        std::cerr << "Warning: n and u direction not orthogonal. Altering u to be consistent with n." << std::endl;
    }

    u = r.cross(n);
    Eigen::Matrix3d ret;
    ret << r[0], u[0], -n[0],
        r[1], u[1], -n[1],
        r[2], u[2], -n[2];
    return ret;
}

int main(int argc, char *argv[])
{
    options opt;
	if(argc < 3 || !parseFlags(argc, argv, opt))
	{
        printUsage();
		return -1;
	}

    m = strtol(argv[2], NULL, 10);

    if (m < 2)
    {
        printUsage();
        return -1;
    }

    if (opt.inflateFactor == std::numeric_limits<double>::infinity())
        opt.inflateFactor = 0.5 / m;

    Eigen::Matrix3d camera = cameraMatrix(opt);

    generateGridPoints(opt.cubeWidth, m, P);

    if (opt.skipgen)
    {
        std::ifstream ifs(argv[1]);
        if (!ifs)
        {
            std::cerr << "Couldn't read file " << argv[1] << std::endl;
            return -1;
        }

        gridvals.resize(m*m*m);
        for (int i = 0; i < m*m*m; i++)
        {
            ifs >> gridvals[i];
        }
    }
    else
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        igl::read_triangle_mesh(argv[1], V, F);
        if(V.rows() == 0)
        {
            std::cerr << "Couldn't read mesh file " << argv[1] << std::endl;
            return -1;
        }

        int nverts = V.rows();
        for (int i = 0; i < nverts; i++)
        {
            V.row(i) = camera.transpose()*V.row(i).transpose();
        }

        meshToSDF(V, F, P, opt.cubeWidth, opt.inflateFactor, 4 * m, gridvals);
    }

    if(opt.outputsdf)
    { 
        std::cout << gridvals << std::endl;
    }
    
    isovalue = 0;

    if (opt.showviewer)
    {
        viewer = new igl::viewer::Viewer;

        recomputeIsosurface();

        //double cubevol = double(m*m*m);

        //std::cout << sdfGradientResidual(opt.cubeWidth, m, gridvals) / cubevol << std::endl;

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
}
