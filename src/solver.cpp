#include "../include/solver.hpp"
#include "../include/types.hpp"
#include "../include/data.hpp"
#include "../include/domain.hpp"
#include <iostream>
#include "../include/stiffnessMatrix.hpp"
//#include <Eigen/PardisoSupport>


using namespace std;

Solver::Solver()
{
 // m_p_data = p_data;
}



void Solver::pcpg(Data& data)
{
//      clock_t begin = clock();

    double eps_iter = 1.0e-4; //atof(options2["eps_iter"].c_str());
    double max_iter = 50; //atoi(options2["max_iter"].c_str());

    //double gPz, gPz_prev, wFw, rho, gamma, norm_gPz0;
    double gPz_prev, wFw, rho, gamma, norm_gPz0;

    Eigen::MatrixXd gtPz, gtPz_prev;

    Eigen::MatrixXd g0, d_rhs, e_loc, e_glb;
    Eigen::MatrixXd iGTG_e, lambda, z, Pz;
    Eigen::MatrixXd Fw, Pg, g, w, w_prev;
    Eigen::MatrixXd beta, alpha;
    Eigen::MatrixXd wtFw;

    Eigen::MatrixXd xx;
    Eigen::MatrixXd yy;
    Eigen::MatrixXd lambda0,invGtG_e;

//    int nSubClst = cluster.get_nSubClst();
//    xx.resize(nSubClst);
//    yy.resize(nSubClst);




    auto &domain = *(data.GetDomain());
    auto &K = *(domain.GetStiffnessMatrix());
    auto rhs_primal = K.GetRHS();

    K.solve(rhs_primal,xx);

    auto  p_opB = data.GetInterfaceOperatorB(); //mult(xx,d_rhs);
    auto  &R_kerK = *(K.GetKernel());

    p_opB->multB(xx, d_rhs);




    // TEST ===============================================
    //Eigen::MatrixXd v1(domain.GetNumberOfPrimalDOFs(),3);
    //v1.setOnes();
    //Eigen::MatrixXd jumps;

    //p_opB->multB(v1, jumps);
    //std::cout <<"jumps: " <<  jumps << '\n';
    // TEST ===============================================


    // !!! e_loc = -Rt*f as well as G = -Rt * Bt
    e_loc = (-1) * R_kerK.transpose() * rhs_primal;

    std::cout << "e_loc = \n" << e_loc << '\n';

    p_opB->mult_invGtG(e_loc,invGtG_e);
    Eigen::MatrixXd RinvGtG_e = R_kerK * invGtG_e;
    p_opB->multB(RinvGtG_e,lambda0);
    lambda0 *= -1;

//    std::cout << "lambda0 = \n" << lambda0 << '\n';


    Eigen::MatrixXd Bt_lambda0;
    p_opB->multBt(lambda0,Bt_lambda0);
    Eigen::MatrixXd Gt_lambda0 = (-1) * R_kerK.transpose() * Bt_lambda0;

    std::cout << "Gt_lambda0 = \n" << Gt_lambda0 << '\n'; 

    Eigen::MatrixXd del = Gt_lambda0 - e_loc;
    std::cout << "|| Gt*lambda0 - e ||  = \n" << del.norm() << '\n';

    MatrixXd Plambda0;
    p_opB->Projection(lambda0,Plambda0);
    std::cout << "Plambda0 = \n" << Plambda0 << '\n';


    // lambda0 = G * inv(GtG) * e
//    iGTG_e.mat_mult_dense(cluster.invGfTGf,"N",e,"N");
//    cluster.mult_Gf(iGTG_e, lambda);

    lambda = lambda0;

    // g0 = F * lambda0 - d_rhs
    p_opB->mult_F(lambda0,g0);
    g0 -= d_rhs;

    g = g0;

    // Pg0
    p_opB->Projection(g,Pg);
    z = Pg;  //    cluster.Preconditioning(Pg,z);
    Pz = z;  //    cluster.Projection(z,Pz,beta);

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
    {
        printf("%4d\t%3.9e \n",it + 1, sqrt(gtPz(0,0)) / norm_gPz0);
    }
        if (sqrt(gtPz(0,0)) < eps_iter * norm_gPz0)
            break;

        p_opB->mult_F(w,Fw);
        wtFw = w.transpose() * Fw;
        domain.hmpi.GlobalSum(wtFw.data(),wtFw.size());
        wtFw *= 0.5;

        rho = -gtPz(0,0) / wtFw(0,0);

        lambda += w * rho;
        g += Fw * rho;

        p_opB->Projection(g,Pg);
        z = Pg;//cluster.Preconditioning(Pg,z);
        Pz = z;//cluster.Projection(z,Pz,beta);

        gtPz_prev = gtPz;
        gtPz = g.transpose() * Pz;
        domain.hmpi.GlobalSum(gtPz.data(),gtPz.size());
        gtPz *= 0.5;

        gamma = gtPz(0,0) / gtPz_prev(0,0);

        w_prev = w;
        w = Pz;
        w = w_prev * gamma;

//        if (options2["vtkWithinIter"].compare("true") == 0)
//            cluster.printVTK(yy, xx, lambda, alpha, it);
    }

#if 0
//    clock_t end = clock();
//    time_solver = double(end - begin) / CLOCKS_PER_SEC;


    // final solution (last parameter -1 avoids numbering of vtk file)
//    printf("Solver time:  %3.1f s.\n",time_solver);
//    printf("Total time:   %3.1f s.\n",time_total);

 //   cluster.printVTK(yy, xx, lambda, alpha, -1);
#endif

}


//Solver::Solver(Mesh* mesh,Data* data)
//{
//  std::vector<int> dirDOFs = mesh->getDirDOFs();
//
//  //  for (auto& id : dirDOFs) cout << id << endl;
//
//  SpMat& K = data->getK();
//  //data->printMatrix(K,"K0.txt"); 
//  data->setDirichlet(K, dirDOFs);
//  //data->printMatrix(K,"K1.txt"); 
//
//
//  Eigen::VectorXd& rhs = data->getRHS();
//
//
//  double sum(0);
//  for (int dof = 0; dof < rhs.size(); dof++)
//    sum += rhs(dof);
//  cout << "sum(rhs) = " << sum - 9.81 << endl;
//
//
//  data->setDirichlet(rhs, dirDOFs);
//
//  m_directSolver.sym_factor(K);
//  m_directSolver.num_factor();
//
//
//  Eigen::VectorXd solution(rhs.size());
//
//  m_directSolver.solve(rhs,solution);
//  mesh->addSolution(solution);
//}

