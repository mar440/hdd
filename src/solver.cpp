#include "../include/solver.hpp"
#include "../include/types.hpp"
#include "../include/data.hpp"
#include "../include/domain.hpp"
#include <iostream>
#include "../include/stiffnessMatrix.hpp"


using namespace std;


bool Solver::pcpg(Data& data, Eigen::VectorXd& _solution)
{


  const int nRHS = 1;

  double eps_iter = 1.0e-5; //atof(options2["eps_iter"].c_str());
  double max_iter = 250; //atoi(options2["max_iter"].c_str());

  double  norm_gPz0;
  
  Eigen::MatrixXd solution, solutionPrew, solutionImK, solutionKerK;
  Eigen::MatrixXd KpBtw, KufBtlambda;

  Eigen::MatrixXd rho(nRHS,nRHS), gamma(nRHS, nRHS);
  Eigen::MatrixXd gtPz, gtPz_prev;

  Eigen::MatrixXd lambda, z, Pz;
  Eigen::MatrixXd Fw, Pg, g, w, w_prev;
  Eigen::MatrixXd alpha, beta, e_loc;
  Eigen::MatrixXd wtFw;

  Eigen::MatrixXd tmp0, tmp;
  Eigen::MatrixXd invGtG_e;

  auto &domain = *(data.GetDomain());
  auto &K = *(domain.GetStiffnessMatrix());
  auto rhs_primal = K.GetRHS();

  auto  p_opB = data.GetInterfaceOperatorB();
  auto  &R_kerK = *(K.GetKernel());

  // !!! e_loc = -Rt*f as well as G = -Rt * Bt
  e_loc = (-1) * R_kerK.transpose() * rhs_primal;

  p_opB->mult_invGtG(e_loc,invGtG_e);
  Eigen::MatrixXd RinvGtG_e = R_kerK * invGtG_e;
  p_opB->multB(RinvGtG_e,lambda);
  lambda *= -1;

  Eigen::MatrixXd Btw;
  p_opB->multBt(lambda,Btw);

#if DBG > 3
  Eigen::MatrixXd Gt_lambda0 = (-1) * R_kerK.transpose() * Btw;
  std::cout << "Gt_lambda0 = \n" << Gt_lambda0 << '\n'; 
  Eigen::MatrixXd del = Gt_lambda0 - e_loc;
  std::cout << "|| Gt*lambda0 - e ||  = \n" << del.norm() << '\n';
  MatrixXd Plambda0;
  beta = p_opB->Projection(lambda,Plambda0);
  std::cout << "Plambda0 = \n" << Plambda0 << '\n';
#endif


  // initial gradient (residual)
  MatrixXd fmBtl = rhs_primal - Btw;
  K.solve(fmBtl,solutionImK);
  p_opB->multB(solutionImK, g);
  g *= -1;

  // solution
  alpha = p_opB->Projection(g,Pg);
  solutionKerK = R_kerK * alpha;
  solution = solutionImK + solutionKerK;

  // equilibrium equation
  K.Mult(solution,KufBtlambda);
  KufBtlambda -= fmBtl;

  // preconditioning
  p_opB->Scaling(Pg);
  p_opB->multBt(Pg,tmp0);
  K.Precond(tmp0,tmp);
  p_opB->multB(tmp,z);
  p_opB->Scaling(z);
  beta = p_opB->Projection(z,Pz);

  gtPz = g.transpose()*Pz;
  domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
  //TODO correction since lambda lives on two subdomains 
  gtPz *= 0.5;

  Eigen::MatrixXd g0tPg0 =  g.transpose() * Pg;
  domain.hmpi.GlobalSum(g0tPg0.data(),g0tPg0.size());
  g0tPg0 *= 0.5;

  double norm_g0Pg0 = sqrt(g0tPg0(0,0));

  norm_gPz0 = sqrt(gtPz(0,0));
  w = Pz;

  if (domain.GetRank() == 0)
  {
    std::cout << "|g0Pg0| = " << norm_g0Pg0 << '\n';
    std::cout << "|gPz0 | = " << norm_gPz0 << '\n';
  }

  if (domain.GetRank() == 0)
  {
    std::cout << "====================================\n";
    //std::cout << " it.   |g(i)|/|g(0)|       |g(i)|\n";
    std::cout << " it. |g(i)|/|g(0)|     K*u-f+Bt*l \n";
    std::cout << "====================================\n";
  }

  int iter(0);
  for (iter = 0; iter < max_iter; iter++)
  {

    if (domain.GetRank() == 0)
    {
      std::cout << std::setw(4) <<  iter+ 1 << "  ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        sqrt(gtPz(0,0)) / norm_gPz0 << "       ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        KufBtlambda.norm() << '\n';
    }


    if (sqrt(gtPz(0,0)) < eps_iter * norm_gPz0 || norm_gPz0 < 1e-14)
      break;



    // F * w
    p_opB->multBt(w, Btw);
    K.solve(Btw,KpBtw);
    p_opB->multB(KpBtw, Fw);

    wtFw = w.transpose() * Fw;
    domain.hmpi.GlobalSum(wtFw.data(),wtFw.size());
    wtFw *= 0.5;

    rho(0,0) = -gtPz(0,0) / wtFw(0,0);

    /////////////
    //  UPDATE
    /////////////

    // dual varibles
    lambda += w * rho;

    // gradient
    g += Fw * rho;

    // primal RHS update
    fmBtl -= Btw * rho;

    // primal variable: Im(K)
    solutionImK -= KpBtw * rho;

    // gradient projection
    alpha = p_opB->Projection(g,Pg);

    // primal varibles: Ker(K)
    solutionKerK = R_kerK * alpha;

    // primal variables
    solution = solutionImK + solutionKerK;

    // primal equlibrium equation
    K.Mult(solution,KufBtlambda);
    KufBtlambda -= fmBtl;

    // preconditioning
    p_opB->Scaling(Pg);
    p_opB->multBt(Pg,tmp0);
    K.Precond(tmp0,tmp);
    p_opB->multB(tmp,z);
    p_opB->Scaling(z);
    beta = p_opB->Projection(z,Pz);

    gtPz_prev = gtPz;
    gtPz = g.transpose() * Pz;
    domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
    gtPz *= 0.5;

    gamma(0,0) = gtPz(0,0) / gtPz_prev(0,0);

    w_prev = w;
    w = Pz;
    w += w_prev * gamma(0,0);

  }



  int niter =  iter + 1;

  std::cout << " HDD solver - finish:\n";
  std::cout << "     number of iteration: " << niter << "\n";

  // copy solution to inout array
  _solution = solution;

  return true;

}

bool Solver::mpcpg(Data& data, Eigen::VectorXd& solution)
{


  auto &domain = *(data.GetDomain());

  int nSubdomains = domain.GetNumberOfSubdomains();

  const int nRHS = nSubdomains;

  double eps_iter = 1.0e-5; //atof(options2["eps_iter"].c_str());
  double max_iter = 250; //atoi(options2["max_iter"].c_str());

  //double gPz, gamma, norm_gPz0;
  //
  double  norm_gPz0;
  Eigen::MatrixXd rho(nRHS,nRHS), gamma(nRHS, nRHS);
  Eigen::MatrixXd gtPz, gtPz_prev;

  Eigen::MatrixXd d_rhs, e_loc;
  Eigen::MatrixXd lambda, z, Pz;
  Eigen::MatrixXd Fw, Pg, g0, g, w, w_prev;
  Eigen::MatrixXd alpha, beta;
  Eigen::MatrixXd wtFw;

  Eigen::MatrixXd xx;
  Eigen::MatrixXd tmp;
  Eigen::MatrixXd invGtG_e;

  // linear opeartor (Kplus)
  auto &K = *(domain.GetStiffnessMatrix());
  // nodal forces - primal RHS
  auto rhs_primal = K.GetRHS();
  // interface operator (B, Bt, G, Gt, inv(GtG), ...)
  auto p_opB = data.GetInterfaceOperatorB();
  // null space of K; R = ker(K)
  auto &R_kerK = *(K.GetKernel());


  // compute initial gradient ( (-1) * residual )
  K.solve(rhs_primal,xx);
  p_opB->multB(xx, d_rhs);
  e_loc = (-1) * R_kerK.transpose() * rhs_primal;

  p_opB->mult_invGtG(e_loc,invGtG_e);
  Eigen::MatrixXd RinvGtG_e = R_kerK * invGtG_e;
  p_opB->multB(RinvGtG_e,lambda);
  lambda *= -1;



  Eigen::MatrixXd Btlambda0;
  p_opB->multBt(lambda,Btlambda0);
  Eigen::MatrixXd Gt_lambda0 = (-1) * R_kerK.transpose() * Btlambda0;

#if DBG > 3
  std::cout << "Gt_lambda0 = \n" << Gt_lambda0 << '\n'; 
  Eigen::MatrixXd del = Gt_lambda0 - e_loc;
  std::cout << "|| Gt*lambda0 - e ||  = \n" << del.norm() << '\n';
  MatrixXd Plambda0;
  beta p_opB->Projection(lambda,Plambda0);
  std::cout << "Plambda0 = \n" << Plambda0 << '\n';
#endif


  // initial gradient (residual)
  // g0 = F * lambda - d_rhs
  p_opB->mult_F(lambda,g0);
  g0 -= d_rhs;
  g = g0;


  beta = p_opB->Projection(g,Pg);
  p_opB->Scaling(Pg);

  p_opB->multBt(Pg,xx);
  K.Precond(xx,tmp);
  p_opB->multB(tmp,z);

  p_opB->Scaling(z);
  beta = p_opB->Projection(z,Pz);


  gtPz = g.transpose()*Pz;
  domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
  //TODO correction since lambda lives on two subdomains 
  gtPz *= 0.5;

  Eigen::MatrixXd g0tPg0 =  g.transpose() * Pg;
  domain.hmpi.GlobalSum(g0tPg0.data(),g0tPg0.size());
  g0tPg0 *= 0.5;

  double norm_g0Pg0 = sqrt(g0tPg0(0,0));

  norm_gPz0 = sqrt(gtPz(0,0));
  w = Pz;

  if (domain.GetRank() == 0)
  {
    std::cout << "|g0Pg0| = " << norm_g0Pg0 << '\n';
    std::cout << "|gPz0 | = " << norm_gPz0 << '\n';
  }

  if (domain.GetRank() == 0)
  {
    std::cout << "====================================\n";
    std::cout << " it.   |g(i)|/|g(0)|       |g(i)|\n";
    std::cout << "====================================\n";
  }

  int iter(0);
  for (iter = 0; iter < max_iter; iter++){

    if (domain.GetRank() == 0)
    {
      std::cout << std::setw(4) <<  iter+ 1 << "  ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(8) <<
        sqrt(gtPz(0,0)) / norm_gPz0 << "  ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(8) <<
        sqrt(gtPz(0,0)) << '\n';
    }


    if (sqrt(gtPz(0,0)) < eps_iter * norm_gPz0 || norm_gPz0 < 1e-14)
      break;

    p_opB->mult_F(w,Fw);
    wtFw = w.transpose() * Fw;
    domain.hmpi.GlobalSum(wtFw.data(),wtFw.size());
    wtFw *= 0.5;

    rho(0,0) = -gtPz(0,0) / wtFw(0,0);

    lambda += w * rho;
    g += Fw * rho;

    beta = p_opB->Projection(g,Pg);
    p_opB->Scaling(Pg);

    p_opB->multBt(Pg,xx);
    K.Precond(xx,tmp);
    p_opB->multB(tmp,z);

    p_opB->Scaling(z);
    beta = p_opB->Projection(z,Pz);

    gtPz_prev = gtPz;
    gtPz = g.transpose() * Pz;
    domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
    gtPz *= 0.5;

    gamma(0,0) = gtPz(0,0) / gtPz_prev(0,0);

    w_prev = w;
    w = Pz;
    w += w_prev * gamma(0,0);

//        if (options2["vtkWithinIter"].compare("true") == 0)
//            cluster.printVTK(yy, xx, lambda, alpha, it);
  }


  {
    Eigen::MatrixXd  Btl, Fl, dFl;
    p_opB->multBt(lambda,Btl);
    Eigen::MatrixXd fmBtl = rhs_primal - Btl;
    Eigen::MatrixXd uI;
    K.solve(fmBtl,uI);
    // g0 = F * lambda - d_rhs
    p_opB->mult_F(lambda,Fl);
    dFl = d_rhs - Fl;
    Eigen::MatrixXd BtdFl, GtdFl;
    p_opB->multBt(dFl,BtdFl);
    GtdFl= (-1) * R_kerK.transpose() * BtdFl;
    p_opB->mult_invGtG(GtdFl,alpha);
    std::cout << "alpha = \n" << alpha << '\n';
    Eigen::MatrixXd uII = R_kerK * alpha;
    solution = uI.col(0) + uII.col(0);
  }

  int niter =  iter + 1;
  
  std::cout << " HDD solver:\n";
  std::cout << " number of iteration: " << niter << "\n";

  return true;

}


