#pragma once

#include <vector>
#include "../include/types.hpp"

namespace linalg {

  void unique(std::vector<int>&);
  Eigen::VectorXd svd0(Eigen::MatrixXd& dmat);
  Eigen::VectorXd svd0(SpMat& smat);
}

namespace tools {

  void printMatrix(SpMat& mat,std::string fname);
  std::vector<int> intersection(
      std::vector<int> v1, std::vector<int> v2);
}



