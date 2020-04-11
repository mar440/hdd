#include "../include/linearAlgebra.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>


namespace linalg {

  void unique(std::vector <int>& vec)
  {
    std::sort(vec.begin(), vec.end());
    std::vector<int>::iterator it;
    it = std::unique(vec.begin(), vec.end());
    vec.resize(distance(vec.begin(), it));
  }


  Eigen::VectorXd svd0(Eigen::MatrixXd& dmat)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dmat,
        Eigen::ComputeThinU | Eigen::ComputeThinV);
    //  cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    //  cout << "Its left singular vectors are the columns of the thin U matrix:" << 
    //    endl << svd.matrixU() << endl;
    //  cout << "Its right singular vectors are the columns of the thin V matrix:" << 
    //    endl << svd.matrixV() << endl;
    //  Vector3f rhs(1, 0, 0);
    //  cout << "Now consider this rhs vector:" << endl << rhs << endl;
    //  cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
    Eigen::MatrixXd singVals = svd.singularValues();
    return singVals;
  }

  Eigen::VectorXd svd0(SpMat& smat)
  {
    Eigen::MatrixXd dmat = Eigen::MatrixXd(smat);
    return svd0(dmat);
  }


  void svd0(const Eigen::MatrixXd& dmat, 
      Eigen::MatrixXd& U, Eigen::VectorXd& diagS, Eigen::MatrixXd& V)
  {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dmat,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    //  cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    //  cout << "Its left singular vectors are the columns of the thin U matrix:" << 
    //    endl << svd.matrixU() << endl;
    //  cout << "Its right singular vectors are the columns of the thin V matrix:" << 
    //    endl << svd.matrixV() << endl;
    //  Vector3f rhs(1, 0, 0);
    //  cout << "Now consider this rhs vector:" << endl << rhs << endl;
    //  cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
   
    U = svd.matrixU();
    diagS = svd.singularValues();
    V = svd.matrixV();
  }



}

namespace tools 
{

  void printMatrix(const SpMat& mat,std::string fname){
    {
      std::ofstream myfile(fname);
      if (myfile.is_open())
      {
        myfile << "\%\%MATRIX MARKET\%\% \n";
        myfile << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() << "\n";
        myfile << std::setprecision(16);
        for (int indJ = 0; indJ < mat.rows(); indJ++) {
          for (SpMat::InnerIterator it(mat, indJ); it; ++it) {
            myfile << it.row() + 1 << " " << it.col() + 1 << " " << it.value() << "\n";
          }
        }
        myfile.close();
      }
    }

  }


  std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2)
  {

    std::vector<int>::iterator it;
    std::vector<int> v12(std::max(v1.size(), v2.size()));

    // std :: set_intersection 
    it= std::set_intersection(v1.begin(), v1.end(),
        v2.begin(), v2.end(), v12.begin()); 
  
    v12.resize(it - v12.begin());
#if DBG > 2
    std::cout << "intersection DOFs (global numbering) " <<std::endl;
    for (auto& iii : v12)
      std::cout << iii << ' ';
    std::cout << std::endl;
#endif
    //std::cout << "The intersection has " << (ls - v12.begin()) << " elements:"; 
    //for (it = v12.begin(); it != ls; ++it) 
    //    std::cout << ' ' << *it; 
    //std::cout << "\n"; 




    return v12;
  }

}
