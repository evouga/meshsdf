#include "SDFAlgs.h"
#include <map>
#include <vector>
#include <set>
#include <deque>
#include <iostream>
#include <limits>
#include <Eigen/Geometry>
#include "igl/ray_mesh_intersect.h"
#include <igl/signed_distance.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/writeOBJ.h>

using namespace std;

struct faceNode
{
    int prevface;
    int nextface;
};

double surfaceArea(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    double totarea = 0;
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d e1 = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d e2 = V.row(F(i, 2)) - V.row(F(i, 0));
        double area = e1.cross(e2).norm();
        totarea += area;
    }
    return totarea;
}

int isInverted(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    if (F.rows() == 0)
        return false;
    // find point guaranteed outside mesh
    double minx = numeric_limits<double>::infinity();
    int nverts = (int)V.rows();
    //igl::writeOBJ("isect.obj", V, F);
    for (int i = 0; i < nverts; i++)
        minx = std::min(minx, V(i, 0));
    Eigen::Vector3d origin(minx - 1.0, 0, 0);

    int nrows = F.rows();
    vector<bool> tried(nrows);
    for (int i = 0; i < nrows; i++)
        tried[i] = false;

    for (int face = 0; face < nrows; face++)
    {
        if (tried[face])
            continue;
        tried[face] = true;

        int curface = face;
        while (true)
        {
            Eigen::Vector3d centroid(0, 0, 0);
            for (int i = 0; i < 3; i++)
            {
                centroid += V.row(F(curface, i));
            }
            centroid /= 3.0;
            Eigen::Vector3d dir = centroid - origin;
            igl::Hit hit;
            bool isect = igl::ray_mesh_intersect(origin, dir, V, F, hit);
            if (!isect)
            {                
                break;
            }
            if (hit.id == curface)
            {
                Eigen::Vector3d e1 = V.row(F(curface, 1)) - V.row(F(curface, 0));
                Eigen::Vector3d e2 = V.row(F(curface, 2)) - V.row(F(curface, 0));
                Eigen::Vector3d n = e1.cross(e2);
                return (n.dot(dir) > 0) ? 1 : 0;
            }
            curface = hit.id;
            tried[curface] = true;
        }
    }
    return -1;
}

void removeInvertedComponents(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF)
{
    map<pair<int, int>, vector<int> > edges;

    int numfaces = (int)F.rows();
    for (int i = 0; i < numfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int curj = F(i,j);
            int nextj = F(i,(j + 1) % 3);
            if (curj > nextj)
                std::swap(curj, nextj);
            edges[pair<int, int>(curj, nextj)].push_back(i);
        }
    }

    bool *visited = new bool[numfaces];
    for (int i = 0; i < numfaces; i++)
        visited[i] = false;

    int componentnum = 0;

    vector<int> keep;

    for (int face = 0; face < numfaces; face++)
    {
        if (visited[face])
            continue;

        //check one connected component
        set<int> component;
        deque<int> q;
        visited[face] = true;
        q.push_front(face);

        while (!q.empty())
        {
            int next = q.front();
            q.pop_front();

            component.insert(next);

            for (int i = 0; i < 3; i++)
            {
                int curi = F(next, i);
                int nexti = F(next, (i + 1) % 3);
                if (curi > nexti)
                    std::swap(curi, nexti);
                for (vector<int>::iterator it = edges[pair<int, int>(curi, nexti)].begin(); it != edges[pair<int, int>(curi, nexti)].end(); ++it)
                {
                    if (visited[*it])
                        continue;
                    q.push_front(*it);
                    visited[*it] = true;
                }
            }
        }

        int componentsize = (int)component.size();
        Eigen::MatrixXi compF(componentsize, 3);
        int idx = 0;
        for (set<int>::iterator it = component.begin(); it != component.end(); ++it)
        {
            compF.row(idx) = F.row(*it);
            idx++;
        }

        double surfacearea = surfaceArea(V, compF);


        if (isInverted(V, compF) == 0)
        {
            // keep entire component
            for (set<int>::iterator it = component.begin(); it != component.end(); ++it)
            {
                keep.push_back(*it);
            }            
        }        

        componentnum++;
    }

    int newsize = (int)keep.size();
    newF.resize(newsize, 3);
    for (int i = 0; i < newsize; i++)
    {
        newF.row(i) = F.row(keep[i]);
    }
}

void reorientFaces(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF)
{
    newF = F;
    map<pair<int, int>, vector<int> > edges;

    int numfaces = (int)F.rows();
    for (int i = 0; i < numfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int curj = F(i,j);
            int nextj = F(i,(j + 1) % 3);
            if (curj > nextj)
                std::swap(curj, nextj);
            edges[pair<int, int>(curj, nextj)].push_back(i);
        }
    }

    bool *visited = new bool[numfaces];
    for (int i = 0; i < numfaces; i++)
        visited[i] = false;

    int componentnum = 0;

    for (int face = 0; face < numfaces; face++)
    {
        if (visited[face])
            continue;

        //orient one connected component
        set<int> component;
        deque<faceNode> q;
        faceNode start;
        start.prevface = -1;
        start.nextface = face;
        visited[face] = true;
        q.push_front(start);

        int flipcount = 0;

        while (!q.empty())
        {
            faceNode next = q.front();
            q.pop_front();

            component.insert(next.nextface);

            if (next.prevface != -1)
            {
                bool flip = false;
                for (int i = 0; i < 3; i++)
                {
                    int nexti = (i + 1) % 3;
                    for (int j = 0; j < 3; j++)
                    {
                        int nextj = (j + 1) % 3;
                        if (newF(next.prevface, i) == newF(next.nextface, j) && newF(next.prevface, nexti) == newF(next.nextface, nextj))
                            flip = true;
                    }
                }
                if (flip)
                {
                    flipcount++;
                    std::swap(newF(next.nextface, 0), newF(next.nextface, 1));
                }
            }

            for (int i = 0; i < 3; i++)
            {
                int curi = F(next.nextface, i);
                int nexti = F(next.nextface, (i + 1) % 3);
                if (curi > nexti)
                    std::swap(curi, nexti);
                for (vector<int>::iterator it = edges[pair<int, int>(curi, nexti)].begin(); it != edges[pair<int, int>(curi, nexti)].end(); ++it)
                {
                    if (visited[*it])
                        continue;
                    faceNode node;
                    node.prevface = next.nextface;
                    node.nextface = *it;
                    q.push_front(node);
                    visited[*it] = true;
                }
            }
        }

        int componentsize = (int)component.size();
        Eigen::MatrixXi compF(componentsize, 3);
        int idx = 0;
        for (set<int>::iterator it = component.begin(); it != component.end(); ++it)
        {
            compF.row(idx) = newF.row(*it);
            idx++;
        }

        if (isInverted(V, compF))
        {
            // reverse entire component
            for (set<int>::iterator it = component.begin(); it != component.end(); ++it)
            {
                std::swap(newF(*it, 0), newF(*it, 1));
            }
            flipcount = componentsize - flipcount;
        }        

        //std::cout << "Flipped " << flipcount << "/" << component.size() << " triangles in component " << componentnum << std::endl;
        
        componentnum++;
    }
}

struct face
{
    int verts[3];

    bool operator<(const face &other) const
    {
        for (int i = 0; i < 3; i++)
        {
            if (verts[i] < other.verts[i])
                return true;
            if (verts[i] > other.verts[i])
                return false;
        }
        return false;
    }
};

void removeDuplicateFaces(const Eigen::MatrixXi &F, Eigen::MatrixXi &newF)
{
    set<face> seen;
    vector<int> keep;
    int nfaces = (int)F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        face f;
        for(int j=0; j<3; j++)
            f.verts[j] = F(i, j);
        std::sort(f.verts, f.verts + 3);
        if (seen.count(f) > 0)
            continue;
        seen.insert(f);
        keep.push_back(i);
    }

    int nnewfaces = (int)keep.size();
    newF.resize(nnewfaces, 3);
    for (int i = 0; i < nnewfaces; i++)
    {
        newF.row(i) = F.row(keep[i]);
    }
}

static int valIndex(int i, int j, int k, int m)
{
    return i + j*m + k*m*m;
}

static double cellGradientNormSquared(int x, int y, int z, double cubewidth, int m, const Eigen::VectorXd &gridvals)
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

void fitIntoCube(Eigen::MatrixXd &V, double cubeWidth)
{
    double mins[3], maxs[3];
    for(int i=0; i<3; i++)
    {
        mins[i] = std::numeric_limits<double>::infinity();
        maxs[i] = -std::numeric_limits<double>::infinity();
    }
    int nverts = (int)V.rows();
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

    // 10% padding
    diameter *= 1.1;

    double s = cubeWidth/diameter;
    if(diameter == 0.0)
        s = 1.0;
    V *= s;

    double offsets[3];
    for (int i = 0; i < 3; i++)
    {
        offsets[i] = 0.5*(cubeWidth - s * (maxs[i] - mins[i]));
    }
    for (int i = 0; i < nverts; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            V(i, j) += offsets[j];
        }
    }
}

void computeSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, Eigen::VectorXd &gridvals)
{

    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    Eigen::MatrixXd N;
    igl::signed_distance(P, V, F, igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), gridvals, I, C, N);
}

void meshToSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, double cubeWidth, double inflateFactor, int subsize, Eigen::VectorXd &gridvals)
{    
    Eigen::MatrixXd rescaledV = V;
    fitIntoCube(rescaledV, cubeWidth);
    Eigen::MatrixXd subP;
    generateGridPoints(cubeWidth, subsize, subP);
    Eigen::VectorXd subgridvals;

    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    Eigen::MatrixXd N;

    double shift = std::max(inflateFactor, 0.5/subsize) * cubeWidth;

    igl::signed_distance(subP, rescaledV, F, igl::SIGNED_DISTANCE_TYPE_UNSIGNED, 0, 3.0*shift, subgridvals, I, C, N);

    for (int i = 0; i < subgridvals.size(); i++)
        if (isnan(subgridvals[i]))
            subgridvals[i] = std::numeric_limits<double>::infinity();

    for (int i = 0; i < subgridvals.size(); i++)
        subgridvals[i] -= shift;
    Eigen::MatrixXd cubeV;
    Eigen::MatrixXi cubeF;
    igl::copyleft::marching_cubes(subgridvals, subP, subsize, subsize, subsize, cubeV, cubeF);

    Eigen::MatrixXi newF;
    removeInvertedComponents(cubeV, cubeF, newF);
    computeSDF(cubeV, newF, P, gridvals);
}