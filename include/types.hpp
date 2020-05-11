#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>




typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> T;
enum PRECONDITIONER_TYPE {NONE, LUMPED, DIRICHLET};

#define GLOBAL_NUMBERING "globalNumbering"
#define SUBDOMAIN_ID "SubdomainId"
#define FORMULATION_ID "FormulationId"
#define MATERIAL_ID "MaterialId"
#define DISPLACEMENT "displacement"

#if DBG > 0
#define HDDTRACES std::cout<< __FILE__ << ':'<<  __LINE__ <<'\n'<< std::flush;
#else
#define HDDTRACES ;
#endif
