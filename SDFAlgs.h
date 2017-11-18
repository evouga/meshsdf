#ifndef SDFALGS_H
#define SDFALGS_H

#include <Eigen/Core>

bool isInverted(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

void removeInvertedComponents(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

void reorientFaces(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

void removeDuplicateFaces(const Eigen::MatrixXi &F, Eigen::MatrixXi &newF);

double sdfGradientResidual(double cubewidth, int m, const Eigen::VectorXd &gridvals);

void generateGridPoints(double cubeWidth, int m, Eigen::MatrixXd &P);

void fitIntoCube(Eigen::MatrixXd &V, double cubeWidth);

void computeSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, Eigen::VectorXd &gridvals);

void meshToSDF(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, double cubeWidth, double inflateFactor, int subsize, Eigen::VectorXd &gridvals);

#endif
