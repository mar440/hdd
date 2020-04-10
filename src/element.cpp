#include "../include/element.hpp"
#include <vtkCell.h>
#include <vtkPoints.h>

using namespace Eigen;
using namespace std;


QUADRATIC_QUAD::QUADRATIC_QUAD() : Element2D() 
{
  m_numberOfGaussPoints = 4;
}

void QUADRATIC_QUAD::shapeFun(VectorXd &N, VectorXd &dNr, VectorXd &dNs, VectorXd &dummyV,
    double r, double s, double dummyD)
{

  N(0) = -0.25*(r - 1.0)*(s - 1.0)*(r + s + 1.0);
  N(1) = -0.25*(r + 1.0)*(s - 1.0)*(r - s - 1.0);
  N(2) = 0.25*(r + 1.0)*(s + 1.0)*(r + s - 1.0);
  N(3) = 0.25*(r - 1.0)*(s + 1.0)*(r - s + 1.0);
  N(4) = 0.5*(r - 1.0)*(r + 1.0)*(s - 1.0);
  N(5) = -0.5*(r + 1.0)*(s - 1.0)*(s + 1.0);
  N(6) = -0.5*(r - 1.0)*(r + 1.0)*(s + 1.0);
  N(7) = 0.5*(r - 1.0)*(s - 1.0)*(s + 1.0);

  dNr(0) = -0.5*r*s + 0.5*r - 0.25*s*s + 0.25*s;
  dNr(1) = -0.5*r*s + 0.5*r + 0.25*s*s - 0.25*s;
  dNr(2) = 0.5*r*s + 0.5*r + 0.25*s*s + 0.25*s;
  dNr(3) = 0.5*r*s + 0.5*r - 0.25*s*s - 0.25*s;
  dNr(4) = 1.0*r*s - 1.0*r;
  dNr(5) = -0.5*s*s + 0.5;
  dNr(6) = -1.0*r*s - 1.0*r;
  dNr(7) = 0.5*s*s - 0.5;

  dNs(0) = -0.25*r*r - 0.5*r*s + 0.25*r + 0.5*s;
  dNs(1) = -0.25*r*r + 0.5*r*s - 0.25*r + 0.5*s;
  dNs(2) = 0.25*r*r + 0.5*r*s + 0.25*r + 0.5*s;
  dNs(3) = 0.25*r*r - 0.5*r*s - 0.25*r + 0.5*s;
  dNs(4) = 0.5*r*r - 0.5;
  dNs(5) = -1.0*r*s - 1.0*s;
  dNs(6) = -0.5*r*r + 0.5;
  dNs(7) = 1.0*r*s - 1.0*s; 
}


void QUADRATIC_QUAD::GaussPoints(int nGP, VectorXd &r, VectorXd &s,
    VectorXd &t, VectorXd &w){

  switch (nGP){
    case 2:
      {
        w << 1.0000000000000000,
        1.0000000000000000 ;
        r << -0.5773502691896257,
          0.5773502691896257 ;
        break;
      }
    case 3:
      {
        w << 0.8888888888888888,
        0.5555555555555556,
        0.5555555555555556 ;
        r << 0.0000000000000000,
          -0.7745966692414834,
          0.7745966692414834 ;
        break;
      }
    case 4:
      {
        w << 0.6521451548625461,
          0.6521451548625461,
          0.3478548451374538,
          0.3478548451374538 ;
        r << -0.3399810435848563,
          0.3399810435848563,
          -0.8611363115940526,
          0.8611363115940526 ;
        break;
      }
  }
}



QUAD::QUAD() : Element2D(){
  m_numberOfGaussPoints = 2;
}

void QUAD::shapeFun(VectorXd &N, VectorXd &dNr, VectorXd &dNs, VectorXd &dummyV, 
    double r, double s, double dummyD){
  N(0) = 0.25*(1.0 - r)*(1.0 - s);
  N(1) = 0.25*(1.0 + r)*(1.0 - s);
  N(2) = 0.25*(1.0 + r)*(1.0 + s);
  N(3) = 0.25*(1.0 - r)*(1.0 + s); 

  dNr(0) = 0.25*(-1.0)*(1.0 - s);
  dNr(1) = 0.25*( 1.0)*(1.0 - s);
  dNr(2) = 0.25*( 1.0)*(1.0 + s);
  dNr(3) = 0.25*(-1.0)*(1.0 + s);

  dNs(0) = 0.25*(1.0 - r)*(-1.0);
  dNs(1) = 0.25*(1.0 + r)*(-1.0);
  dNs(2) = 0.25*(1.0 + r)*( 1.0);
  dNs(3) = 0.25*(1.0 - r)*( 1.0);
}


void QUAD::GaussPoints(int nGP, VectorXd &r, VectorXd &s, VectorXd &t, VectorXd &w){

  switch (nGP){
    case 2:
      {
        w << 1.0000000000000000,
        1.0000000000000000 ;
        r << -0.5773502691896257,
          0.5773502691896257 ;
        break;
      }
    case 3:
      {
        w << 0.8888888888888888,
        0.5555555555555556,
        0.5555555555555556 ;
        r << 0.0000000000000000,
          -0.7745966692414834,
          0.7745966692414834 ;
        break;
      }
    case 4:
      {
        w << 0.6521451548625461,
          0.6521451548625461,
          0.3478548451374538,
          0.3478548451374538 ;
        r << -0.3399810435848563,
          0.3399810435848563,
          -0.8611363115940526,
          0.8611363115940526 ;
        break;
      }
  }
//	switch (nGP){
//	case 4:
//	{
//		w << 1.0000000000000000,
//			1.0000000000000000,
//			1.0000000000000000,
//			1.0000000000000000;
//			r << -0.5773502691896257,
//				0.5773502691896257,
//				-0.5773502691896257,
//				0.5773502691896257;
//			s << -0.5773502691896257,
//				-0.5773502691896257,
//				0.5773502691896257,
//				0.5773502691896257;
//		break;
//	}
//	case 9:
//	{
//		w << 0.8888888888888888,
//			0.5555555555555556,
//			0.5555555555555556;
//		r << 0.0000000000000000,
//			-0.7745966692414834,
//			0.7745966692414834;
//		break;
//	}
//	case 16:
//	{
//		w << 0.6521451548625461,
//			0.6521451548625461,
//			0.3478548451374538,
//			0.3478548451374538;
//		r << -0.3399810435848563,
//			0.3399810435848563,
//			-0.8611363115940526,
//			0.8611363115940526;
//		break;
//	}
//	}
}





void Element2D::assembly_elasticity(MatrixXd& K, VectorXd &f, vtkCell* cell, double E_mu[])
{

  int nP = cell->GetNumberOfPoints();
  const int dim = 2;

//  auto CellType = cell->GetCellType(); 
  vtkPoints *nodesOfElem = cell->GetPoints(); 
  Matrix3d C = Matrix3d::Zero();
  MatrixXd B = MatrixXd::Zero(3, nP * dim);
  K.resize(nP * dim, nP * dim);
  K.setZero();
  f.resize(nP * dim);
  f.setZero();
  double E = E_mu[0]; 
  double mu = E_mu[1];
  double rho = E_mu[2];
  C(0, 0) = 1.0;
  C(1, 1) = 1.0;
  C(0, 1) = mu;
  C(1, 0) = mu;
  C(2, 2) = 1.0 - mu; 
  C *=  E / (1.0 - mu * mu); 


  VectorXd N(nP);
  VectorXd dNr(nP);
  VectorXd dNs(nP);
  VectorXd dummyVec;
  VectorXd r_GP(m_numberOfGaussPoints);
  VectorXd s_GP(m_numberOfGaussPoints);
  VectorXd t_GP(m_numberOfGaussPoints);
  VectorXd weights(m_numberOfGaussPoints);
  Vector3d volF;

  volF(0) = 0; 
  volF(1) = -9.81;

  GaussPoints(m_numberOfGaussPoints, r_GP, s_GP, t_GP, weights);
  double x_rs, y_rs, dxdr, dxds, dydr, dyds;



  MatrixXd dNrs(dim, nP);
  MatrixXd dNxy(dim, nP);
  Matrix2d J; 
  double dummyD = 0;
  double _x, _y;



  for (int i = 0; i < m_numberOfGaussPoints; i++){
    for (int j = 0; j < m_numberOfGaussPoints; j++){
      shapeFun(N, dNr, dNs, dummyVec, r_GP(i), r_GP(j), dummyD);
      x_rs = y_rs = dxdr = dxds = dydr = dyds = 0;

      for (int k = 0; k < nP; k++){
        auto *i_xyz = nodesOfElem->GetPoint(k);
        _x = i_xyz[0];
        _y = i_xyz[1];

        x_rs += N(k) * _x;
        y_rs += N(k) * _y;

        dxdr += dNr(k) * _x;
        dydr += dNr(k) * _y;

        dxds += dNs(k) * _x;
        dyds += dNs(k) * _y;
      }



      J(0, 0) = dxdr;
      J(1, 0) = dxds;
      J(0, 1) = dydr;
      J(1, 1) = dyds; 
      double detJ = (J.determinant()); 


      if (detJ < 0)
        cout << "negativ Jacobian" << endl;

      for (int k = 0; k < nP; k++){
        dNrs(0, k) = dNr(k);
        dNrs(1, k) = dNs(k);
      }
      dNxy = J.inverse() * dNrs; 

      for (int k = 0; k < nP; k++){
        B(0, k) = dNxy(0, k);
        B(1, k + nP) = dNxy(1, k);
        B(2, k) = dNxy(1, k);
        B(2, k + nP) = dNxy(0, k);
      }
      double w_ij_detJ = weights(i) * weights(j) * detJ;
      K.noalias() += (B.transpose() * (C * B) * w_ij_detJ);
      for (int iF = 0; iF < dim; iF++)
        f.segment(iF * nP, nP) += N.transpose() * rho *  (volF(iF) * w_ij_detJ); 
    }
  }
}
