#ifndef SDFALGS_H
#define SDFALGS_H

#include <Eigen/Core>

// 0: not inverted
// 1: inverted
// -1: cannot determine (due to size too small, or other problem with ray tracing)
int isInverted(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

void removeInvertedComponents(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

void reorientFaces(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

void removeDuplicateFaces(const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

double sdfGradientResidual(double cubewidth, int m, const Eigen::VectorXd &gridvals);

void generateGridPoints(double cubeWidth, int m, Eigen::MatrixXd &P);

void fitIntoCube(Eigen::MatrixXd &V, double cubeWidth);

void computeSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, Eigen::VectorXd &gridvals);

void meshToSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, double cubeWidth, double inflateFactor, int subsize, Eigen::VectorXd &gridvals);

#endif
