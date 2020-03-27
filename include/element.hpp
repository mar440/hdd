#pragma once 

#include <Eigen/Dense>

class vtkCell;


using namespace Eigen;

class Element{
  public:
    Element(){}
    ~Element(){}
    virtual void assembly_elasticity(MatrixXd& K, VectorXd&, 
        vtkCell* cell, double E_mu[]) = 0;
    virtual void GaussPoints(int nGP, VectorXd &r, VectorXd &s,
        VectorXd &t, VectorXd &w) = 0;
    virtual void shapeFun(VectorXd &N, VectorXd &dNr, VectorXd &dNs,
        VectorXd &dNt, double r, double s, double t) = 0;
  protected:
    int m_numberOfGaussPoints=0;
};



class Element2D : public Element 
{
  private:
    int m_nOfGP;
  public:
    Element2D(){}
    Element2D(int _nOfGP) :  m_nOfGP(_nOfGP) {}
    void assembly_elasticity(MatrixXd &, VectorXd& f, vtkCell*, double[]);
};




class QUADRATIC_QUAD : public Element2D {
  public:
    QUADRATIC_QUAD();
    QUADRATIC_QUAD(int _nGP) : Element2D(_nGP){}
    ~QUADRATIC_QUAD(); 
    void GaussPoints(int nGP, VectorXd &r, VectorXd &s, VectorXd &t, VectorXd &w);
    void shapeFun(VectorXd &N, VectorXd &dNr, VectorXd &dNs, VectorXd &dNt, double r, double s, double t);

};
