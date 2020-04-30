#pragma once

#include <vector>
#include "../include/types.hpp"

namespace linalg {

  void unique(std::vector<int>&);
  Eigen::VectorXd svd0(Eigen::MatrixXd& dmat);
  Eigen::VectorXd svd0(SpMat& smat);

  void svd0(const Eigen::MatrixXd& dmat,
      Eigen::MatrixXd& U,
      Eigen::VectorXd& diagS,
      Eigen::MatrixXd& V);

}

namespace tools {

  void printMatrix(const SpMat& mat,std::string fname);
  void printMatrix(const Eigen::MatrixXd& mat,std::string fname);
  std::vector<int> intersection(
      std::vector<int> v1, std::vector<int> v2);
//  template<typename T>
    void printArray(std::vector<int>& input, std::string fname);
}



