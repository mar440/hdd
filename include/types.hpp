#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>




typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> T;


#define GLOBAL_NUMBERING "globalNumbering"
#define SUBDOMAIN_ID "SubdomainId"
#define FORMULATION_ID "FormulationId"
#define MATERIAL_ID "MaterialId"
#define DISPLACEMENT "displacement"
