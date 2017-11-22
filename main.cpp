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
};

bool parseFlags(int argc, char *argv[], options &opt)
{
    opt.skipgen = false;
    opt.outputsdf = false;
    opt.showviewer = false;
    opt.cubeWidth = 1.0;
    opt.inflateFactor = std::numeric_limits<double>::infinity();
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
