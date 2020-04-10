#include "../include/solver.hpp"
#include "../include/types.hpp"
#include "../include/data.hpp"
#include "../include/domain.hpp"
#include <iostream>
#include "../include/stiffnessMatrix.hpp"


using namespace std;

Solver::Solver()
{
 // m_p_data = p_data;
}



void Solver::pcpg(Data& data, Eigen::VectorXd& solution)
{
//      clock_t begin = clock();

    double eps_iter = 1.0e-5; //atof(options2["eps_iter"].c_str());
    double max_iter = 250; //atoi(options2["max_iter"].c_str());

    //double gPz, wFw, rho, gamma, norm_gPz0;
    double rho, gamma, norm_gPz0;

    Eigen::MatrixXd gtPz, gtPz_prev;

    Eigen::MatrixXd g0, d_rhs, e_loc, e_glb;
    Eigen::MatrixXd iGTG_e, lambda, z, Pz;
    Eigen::MatrixXd Fw, Pg, g, w, w_prev;
    Eigen::MatrixXd beta, alpha;
    Eigen::MatrixXd wtFw;

    Eigen::MatrixXd xx;
    Eigen::MatrixXd tmp;
    Eigen::MatrixXd invGtG_e;





    auto &domain = *(data.GetDomain());
    auto &K = *(domain.GetStiffnessMatrix());
    auto rhs_primal = K.GetRHS();

    K.solve(rhs_primal,xx);

    auto  p_opB = data.GetInterfaceOperatorB(); //mult(xx,d_rhs);
    auto  &R_kerK = *(K.GetKernel());

    p_opB->multB(xx, d_rhs);



    // !!! e_loc = -Rt*f as well as G = -Rt * Bt
    e_loc = (-1) * R_kerK.transpose() * rhs_primal;


    p_opB->mult_invGtG(e_loc,invGtG_e);
    Eigen::MatrixXd RinvGtG_e = R_kerK * invGtG_e;
    p_opB->multB(RinvGtG_e,lambda);
    lambda *= -1;



    Eigen::MatrixXd Bt_lambda0;
    p_opB->multBt(lambda,Bt_lambda0);
    Eigen::MatrixXd Gt_lambda0 = (-1) * R_kerK.transpose() * Bt_lambda0;

#if DBG > 3
    std::cout << "Gt_lambda0 = \n" << Gt_lambda0 << '\n'; 
    Eigen::MatrixXd del = Gt_lambda0 - e_loc;
    std::cout << "|| Gt*lambda0 - e ||  = \n" << del.norm() << '\n';
    MatrixXd Plambda0;
    p_opB->Projection(lambda,Plambda0);
    std::cout << "Plambda0 = \n" << Plambda0 << '\n';
#endif


    // initial gradient (residual)
    // g0 = F * lambda - d_rhs
    p_opB->mult_F(lambda,g0);
    g0 -= d_rhs;
    g = g0;


    p_opB->Projection(g,Pg);
    p_opB->Scaling(Pg);

    p_opB->multBt(Pg,xx);
    K.mult(xx,tmp);
    //tmp = xx;  
    p_opB->multB(tmp,z);

    p_opB->Scaling(z);
    p_opB->Projection(z,Pz);


    gtPz = g.transpose()*Pz;
    domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
    //TODO correction since lambda lives on two subdomains 
    gtPz *= 0.5;

    Eigen::MatrixXd g0tPg0 =  g.transpose() * Pg;
    domain.hmpi.GlobalSum(g0tPg0.data(),g0tPg0.size());
    g0tPg0 *= 0.5;

    double norm_g0Pg0 = sqrt(g0tPg0(0,0));

    norm_gPz0 = sqrt(gtPz(0,0));

    if (domain.GetRank() == 0)
    {
      printf("\n|g0Pg0| = %3.9e  \n", norm_g0Pg0);
      printf(  "|gPz0|  = %3.9e  \n\n", norm_gPz0);
    }
    w = Pz;

    if (domain.GetRank() == 0)
    {
      printf("=======================\n");
      printf(" it.\t||gradient||\n");
      printf("=======================\n");
    }

    for (int it = 0; it < max_iter; it++){

      if (domain.GetRank() == 0)
        printf("%4d\t%3.9e \n",it + 1, sqrt(gtPz(0,0)) / norm_gPz0);

      if (sqrt(gtPz(0,0)) < eps_iter * norm_gPz0 || norm_gPz0 == 0)
        break;

      p_opB->mult_F(w,Fw);
      wtFw = w.transpose() * Fw;
      domain.hmpi.GlobalSum(wtFw.data(),wtFw.size());
      wtFw *= 0.5;

      rho = -gtPz(0,0) / wtFw(0,0);

      lambda += w * rho;
      g += Fw * rho;

      p_opB->Projection(g,Pg);
      p_opB->Scaling(Pg);

      p_opB->multBt(Pg,xx);
      K.mult(xx,tmp);
      //tmp = xx;  
      p_opB->multB(tmp,z);

      p_opB->Scaling(z);
      p_opB->Projection(z,Pz);

      gtPz_prev = gtPz;
      gtPz = g.transpose() * Pz;
      domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
      gtPz *= 0.5;

      gamma = gtPz(0,0) / gtPz_prev(0,0);

      w_prev = w;
      w = Pz;
      w += w_prev * gamma;

//        if (options2["vtkWithinIter"].compare("true") == 0)
//            cluster.printVTK(yy, xx, lambda, alpha, it);
    }


    Eigen::MatrixXd  Btl, Fl, dFl;
    p_opB->multBt(lambda,Btl);

    Eigen::MatrixXd f_Btl = rhs_primal - Btl;

    Eigen::MatrixXd uI;

    K.solve(f_Btl,uI);


    // g0 = F * lambda - d_rhs
    p_opB->mult_F(lambda,Fl);
    dFl = d_rhs - Fl;

    Eigen::MatrixXd BtdFl, GtdFl;
    p_opB->multBt(dFl,BtdFl);
    GtdFl= (-1) * R_kerK.transpose() * BtdFl;

    p_opB->mult_invGtG(GtdFl,alpha);

    Eigen::MatrixXd uII = R_kerK * alpha;

    solution = uI.col(0) + uII.col(0);




}


