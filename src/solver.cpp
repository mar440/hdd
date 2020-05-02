#include "../include/solver.hpp"
#include "../include/types.hpp"
#include "../include/data.hpp"
#include "../include/domain.hpp"
#include "../include/stiffnessMatrix.hpp"

#include <iostream>
#include <chrono>




using namespace std;


bool Solver::pcpg(Data& data, Eigen::VectorXd& _solution)
{


  auto startTime = std::chrono::steady_clock::now();

  auto solverOpts = data.GetChild("solver");
  int verboseLevel = data.GetChild("outputs").get<int>("verbose");

  double eps_iter = solverOpts.get<double>("stopingCriteria",0);
  int max_iter= solverOpts.get<int>("maxNumbIter",0);


  std::cout <<"eps_iter: " << eps_iter << '\n';
  std::cout <<"max_iter: " << max_iter << '\n';



  const int nRHS = 1;


  double  norm_gPz0;

  // primal dim
  Eigen::MatrixXd solution, solutionImK, solutionKerK;
  Eigen::MatrixXd solutionPrev, del_solution;
  Eigen::MatrixXd KpBtw, KufBtlambda;

  // dual dim
  Eigen::MatrixXd lambda, z, Pz;
  Eigen::MatrixXd Fw, Pg, g, w, w_prev, Btw;

  // trial dim
  Eigen::MatrixXd alpha, beta, e_loc;

  // scalars
  Eigen::MatrixXd rho(nRHS,nRHS), gamma(nRHS, nRHS);
  Eigen::MatrixXd wtFw, gtPz, gtPz_prev;

  // tmps
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

  p_opB->multBt(lambda,Btw);

  if (verboseLevel > 3){
    Eigen::MatrixXd Gt_lambda0 = (-1) * R_kerK.transpose() * Btw;
    std::cout << "Gt_lambda0 = \n" << Gt_lambda0 << '\n'; 
    Eigen::MatrixXd del = Gt_lambda0 - e_loc;
    std::cout << "|| Gt*lambda0 - e ||  = \n" << del.norm() << '\n';
    MatrixXd Plambda0;
    beta = p_opB->Projection(lambda,Plambda0);
    std::cout << "Plambda0 = \n" << Plambda0 << '\n';
  }

  // initial gradient (residual)
  MatrixXd fmBtl = rhs_primal - Btw;
  K.solve(fmBtl,solutionImK);
  p_opB->multB(solutionImK, g);
  g *= -1;

  // solution
  alpha = p_opB->Projection(g,Pg);
  solutionKerK = R_kerK * alpha;
  solution = solutionImK + solutionKerK;

  solutionPrev = Eigen::MatrixXd::Zero(rhs_primal.rows(),nRHS);
  del_solution = solution;
 
  double norm_solution0 = pow(del_solution.norm(),2);

  domain.hmpi.GlobalSum(&norm_solution0,1);
  norm_solution0 = sqrt(norm_solution0);



  double norm_solution(1);

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

    std::cout << "\n\n\n Relative norms: ";
    std::cout << "\n 1- projected precond. dual gradient    P*(F*lam(i)-d),";
    std::cout << "\n 2- primal solution difference          u(i)-u(i-1)";
    std::cout << "\n 2- primal gradient                     K*u-f+Bt*lam(i)";
    std::cout << "\n===========================================================";
    std::cout << "\n it. |g(i)|/|g(0)|    |du| / |u0|     K*u-f+Bt*l";
    std::cout << "\n===========================================================\n";
  }


  int iter(0);
  for (iter = 0; iter < max_iter; iter++)
  {

    norm_solution = pow(del_solution.norm(),2);
    domain.hmpi.GlobalSum(&norm_solution,1);
    norm_solution = sqrt(norm_solution);

    if (domain.GetRank() == 0)
    {
      std::cout << std::setw(4) <<  iter+ 1 << "  ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        sqrt(gtPz(0,0)) / norm_gPz0 << "      ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        norm_solution / norm_solution0 << "      ";
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

    del_solution = solution - solutionPrev;
    solutionPrev = solution;

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




    auto endTime = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime-startTime;
    std::cout << std::fixed << std::setprecision(2) << 
      "Solver time: " << elapsed_seconds.count() << " s\n";



  return true;

}


bool Solver::mpcpg(Data& data, Eigen::VectorXd& _solution)
{


  const int nRHS = 1;

  double eps_iter = 1.0e-5; //atof(options2["eps_iter"].c_str());
  double max_iter = 250; //atoi(options2["max_iter"].c_str());

  double  norm_gPz0;

  // primal dim
  Eigen::MatrixXd solution, solutionImK, solutionKerK;
  Eigen::MatrixXd solutionPrev, del_solution;
  Eigen::MatrixXd KpBtw, KufBtlambda;

  // dual dim
  Eigen::MatrixXd lambda, z, Pz;
  Eigen::MatrixXd Fw, Pg, g, w, w_prev, Btw;

  // trial dim
  Eigen::MatrixXd alpha, beta, e_loc;

  // scalars
  Eigen::MatrixXd rho(nRHS,nRHS), gamma(nRHS, nRHS);
  Eigen::MatrixXd wtFw, gtPz, gtPz_prev;

  // tmps
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

  solutionPrev = Eigen::MatrixXd::Zero(rhs_primal.rows(),nRHS);
  del_solution = solution;
 
  double norm_solution0 = pow(del_solution.norm(),2);

  domain.hmpi.GlobalSum(&norm_solution0,1);
  norm_solution0 = sqrt(norm_solution0);



  double norm_solution(1);

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

// SFETI
  p_opB->SFETI_Beta(Pz);


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

    std::cout << "\n\n\n Relative norms: ";
    std::cout << "\n 1- projected precond. dual gradient    P*(F*lam(i)-d),";
    std::cout << "\n 2- primal solution difference          u(i)-u(i-1)";
    std::cout << "\n 2- primal gradient                     K*u-f+Bt*lam(i)";
    std::cout << "\n===========================================================";
    std::cout << "\n it. |g(i)|/|g(0)|    |du| / |u0|     K*u-f+Bt*l";
    std::cout << "\n===========================================================\n";
  }


  int iter(0);
  for (iter = 0; iter < max_iter; iter++)
  {

    norm_solution = pow(del_solution.norm(),2);
    domain.hmpi.GlobalSum(&norm_solution,1);
    norm_solution = sqrt(norm_solution);

    if (domain.GetRank() == 0)
    {
      std::cout << std::setw(4) <<  iter+ 1 << "  ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        sqrt(gtPz(0,0)) / norm_gPz0 << "      ";
      std::cout << std::setw(4) << std::scientific << std::setprecision(4) <<
        norm_solution / norm_solution0 << "      ";
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

    del_solution = solution - solutionPrev;
    solutionPrev = solution;

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



