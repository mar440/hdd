#include "../include/stiffnessMatrix.hpp"
#include "../include/linearAlgebra.hpp"
#include <math.h>
#include <stdexcept>
#include <algorithm>
#include <ctime>
#include <Eigen/Householder>

StiffnessMatrix::StiffnessMatrix(
    std::map<int,int>* pg2l,int my_rank)
{

  m_pg2l = pg2l;

  m_cnt_setLocalMatrix =0;
  m_cnt_setLocalRHS=0;

  m_numberOfElements = 0;
  m_spmatK.resize(0,0);
  m_trK.resize(0);
  m_myrank = my_rank;
  m_neq = m_pg2l->size();
  m_rhs.resize(0);
  //m_preconditionerType = LUMPED;
  m_preconditionerType = DIRICHLET;

  m_K_IB.resize(0,0);
  m_K_BI.resize(0,0);
  m_K_BB.resize(0,0);

  //m_preconditionerType = NONE;

//    SetDirichletPrecond();


}



void StiffnessMatrix::AddElementContribution(std::vector<int>& glbIds,
    std::vector<double>& valLocK)
{

  //TODO if patern does not change
  // use exixting triplets and replace Value() only

  if (m_cnt_setLocalMatrix==0)
  {
//RHS    m_rhs.resize(m_neqPrimal);
//RHS    m_rhs.setZero();
#if DBG > 0
    std::cout << "local numbering - start\n";
#endif
  }


  int neqLocal = glbIds.size();

  if (pow(neqLocal,2) != valLocK.size())
    std::runtime_error(__FILE__);

  int rowLocal(0);

  for (int row = 0; row < neqLocal; row++)
  {
    rowLocal = (*m_pg2l)[glbIds[row]];
//RHS    m_rhs(rowLocal) += valLocRHS[row];

    for (int col = 0; col < neqLocal; col++)
    {
      m_trK.push_back(
          T(rowLocal, (*m_pg2l)[glbIds[col]],
            valLocK[row + neqLocal * col]));
    }
#if DBG > 3
    std::cout << (*m_pg2l)[glbIds[row]] << ' ';
//  m_pcomm =  _pcomm;
#endif
  }
#if DBG > 3
  std::cout << '\n';
#endif

  m_cnt_setLocalMatrix++;

}


void StiffnessMatrix::AddRHSContribution(std::vector<int>& glbIds,
    std::vector<double>& valLocF)
{

  if (m_cnt_setLocalRHS==0)
  {
    m_rhs.resize(m_neq);
    m_rhs.setZero();
#if DBG > 0
    std::cout << "local numbering (RHS) - start\n";
#endif
  }



  int neqLocal = glbIds.size();

  if (pow(neqLocal,2) != valLocF.size())
    std::runtime_error(__FILE__);

  int rowLocal(0);

  for (int row = 0; row < neqLocal; row++)
  {
    rowLocal = (*m_pg2l)[glbIds[row]];
    m_rhs(rowLocal) += valLocF[row];
  }

  m_cnt_setLocalRHS++;

}

void StiffnessMatrix::m_FactorizeLinearOperator(std::vector<int> nullPivots)
{

  double penalty(0);
  if (nullPivots.size() > 0)
  {
    auto diags = m_spmatK.diagonal();
    for (int iv = 0; iv < diags.size(); iv++)
      penalty += diags(iv);

    penalty /= diags.size();

    for (auto& pvt : nullPivots)
      m_spmatK.coeffRef(pvt,pvt) +=  penalty;

  }

  m_pardisoSolver.analyzePattern(m_spmatK);
  m_pardisoSolver.factorize(m_spmatK);

  if (nullPivots.size() > 0)
  {
    for (auto& pvt : nullPivots)
      m_spmatK.coeffRef(pvt,pvt) -=  penalty;
  }


}



void StiffnessMatrix::FinalizeNumericPart(const std::vector<int>& DirDOFs)
{


  if (m_numberOfElements==0) m_numberOfElements = m_cnt_setLocalMatrix;

  m_spmatK.resize(m_neq, m_neq);
  m_spmatK.setFromTriplets(m_trK.begin(), m_trK.end());
  m_spmatK.makeCompressed();
  ApplyDirichletBC(m_spmatK,DirDOFs);


  m_trK.clear();
  m_trK.shrink_to_fit();



  // build-in factorization based on PARDISO
  m_kerK = _GetKernelFromK(m_spmatK);
  m_nullPivots = GetNullPivots(m_kerK);
  m_FactorizeLinearOperator(m_nullPivots);

#if DBG>2
  Eigen::MatrixXd v2 =  Eigen::MatrixXd::Random(m_neq,2);

  Eigen::MatrixXd Kv2 = m_spmatK * v2;
  Eigen::MatrixXd KpKv2 = m_pardisoSolver.solve(Kv2);
  Eigen::MatrixXd KKpKv2 = m_spmatK * KpKv2;
  Eigen::MatrixXd del = KKpKv2 - Kv2;
  std::cout << "sol =  " << del.norm() << std::endl;
#endif



#if DBG>2
  Eigen::Map<Eigen::VectorXi> MnullPivots(m_nullPivots.data(),m_nullPivots.size());
  std::cout << "nullPivots: \n"  <<
    MnullPivots << std::endl;
#endif


  


#if DBG>3
  m_dbg_printStiffnessMatrix(m_spmatK);
#endif


}


void StiffnessMatrix::Precond(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{
  
  if (in.rows() != m_spmatK.rows())
    std::runtime_error(__FILE__);

  if (m_preconditionerType == NONE)
    out = in;
  else if (m_preconditionerType == LUMPED)
    out = m_spmatK * in;
  else if (m_preconditionerType == DIRICHLET)
  {
    Eigen::MatrixXd inB(m_B_DOFs.size(),in.cols());

    for (int dof = 0; dof < (int) m_B_DOFs.size(); dof++)
      inB.row(dof) = in.row(m_B_DOFs[dof]);

    Eigen::MatrixXd outB = m_K_BB * inB;
    Eigen::MatrixXd rhsB  = m_K_IB * inB;
    outB -= m_K_BI * (m_pardisoSolverDirichlet_KII.solve(rhsB));

    out.resize(in.rows(), in.cols()); out.setZero();

    for (int dof = 0; dof < (int) m_B_DOFs.size(); dof++)
      out.row(m_B_DOFs[dof]) = outB.row(dof);

  }
  else
    std::runtime_error("UNKNOWN PRECONDITIONER TYPE.");


}


void StiffnessMatrix::Mult(const Eigen::MatrixXd& in, Eigen::MatrixXd& out )
{
  if (in.rows() != m_spmatK.rows())
    std::runtime_error(__FILE__);

  out = m_spmatK * in;


}

void StiffnessMatrix::solve(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{
  out = m_pardisoSolver.solve(in);
}




Eigen::MatrixXd StiffnessMatrix::_GetKernelFromK(const SpMat& K)
{

//
// Routine calculates kernel Kplus_R of K satisfied euqality K * Kplus_R = O,
// where O is zero matrix, and it makes the matrix K non-singular (K_reg)
// utilizing spectral conditions of Schur complement. Then ||K-K*inv(K_reg)*K||=0.0
//
//
// rev. 2016-02-03 (AM)
// rev. 2017-06-23 (AM)
// rev. 2020-04-02 (AM)
//==============================================================================
//

typedef int  eslocal;
#define SEQ_VECTOR std::vector
//
//    1) diagonalScaling
//  reducing of big jump coefficient effect (TODO include diagonal scaling into whole ESPRESO)
//BOOL DIAGONALSCALING                                  = true;
  bool diagonalScaling                                  = true;

//    2) permutVectorActive
//  random selection of singular DOFs
// 0 - no permut., 1 - std::vector shuffle, 2 - generating own random sequence -
//ESLOCAL PERMUTVECTORACTIVE                            = 1;
  eslocal permutVectorActive                            = 1;

//    3) use_null_pivots_or_s_set
  // NtN_Mat from null pivots or fixing DOFs
//BOOL USE_NULL_PIVOTS_OR_S_SET                         = TRUE;
  bool use_null_pivots_or_s_set                         = true;

////    4) diagonalRegularization
////  regularization only on diagonal elements (big advantage: patern of K and K_regular is the same !!!)
////  size of set 's' = defect(K)
////  It's is active, only if and only if 'use_null_pivots_or_s_set = true'
////BOOL DIAGONALREGULARIZATION                           = TRUE;
//  bool diagonalRegularization                           = true;


//    6) get_n_first_and_n_last_eigenvals_from_dense_S
// get and print 2*n S eigenvalues
//ESLOCAL GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_S = 0;
  eslocal get_n_first_and_n_last_eigenvals_from_dense_S = 0;


//    8) fixing_nodes_or_dof
// non-singular part determined by fixing nodes (FN),
// min(fixing_nodes_or_dof)>=3; if variable is nonzero,
// parameter sc_size is set to fixing_nodes_or_dof*dofPerNode
//ESLOCAL FIXING_NODES_OR_DOF                           = 0;
  eslocal fixing_nodes_or_dof = 0;
//ESLOCAL DOFPERNODE                                    = 3;
  eslocal dofPerNode                                    = 3;
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//

//    2) check_nonsing
// if check_nonsing>0, checking of K_rr non-singularity is activated and it is repeated
// (check_nonsing) times.
//ESLOCAL CHECK_NONSING                                 = 0;

//    3) max_size_of_dense_matrix_to_get_eigs
// if size of K is less then CHECK_N..., K is converted to dense format to get eigenvalues.
//ESLOCAL MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS          = 2500;

//    4) sc_size
// specification of size of Schur complement used for detection of zero eigenvalues.
//eslocal  sc_size >= expected defect 'd' (e.g. in elasticity d=6).
//ESLOCAL SC_SIZE                                       = 50;
  eslocal sc_size0                                      = 40;

//    5) twenty
// testing last twenty eigenvalues of S to distinguish, if d-last ones are zero or not.
//ESLOCAL TWENTY                                        = 20;
  eslocal twenty                                        = 20;
  // twenty eigenvalues are ascendly ordered in d = d[0],d[1], ..., d[n-2],d[n-1]

//    6) jump_in_eigenvalues_alerting_singularity
// if d[i]/d[i+1]< jump_in_eigenvalues_alerting_singularity, d[i] is last nonzero eigenvalue
//DOUBLE JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY       = 1.0E-5;
  double jump_in_eigenvalues_alerting_singularity       = 1.0e-5;

// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////






  int K_rows = K.rows();

//  if (!use_null_pivots_or_s_set) diagonalRegularization=false;


  //TODO if K.rows<=sc_size, use directly input K instead of S
  //
  int n_nodsSub = 0;
  if (fixing_nodes_or_dof>0){
    sc_size0 = fixing_nodes_or_dof*dofPerNode;
    n_nodsSub = round(K_rows / dofPerNode);
  }
  //
  //##########################################################################################
  //
  eslocal i_start = 0;

  eslocal const sc_size  = K_rows < sc_size0 ? K_rows : sc_size0;

  if (K_rows < sc_size0) fixing_nodes_or_dof=0;


  eslocal const nonsing_size = K_rows - sc_size - i_start;
  SEQ_VECTOR <eslocal > permVec;

#if DBG>2
  eslocal j_start = nonsing_size;
  std::cout << "nonsing_size = " << nonsing_size << std::endl;
  std::cout << "j_start      = " << j_start << std::endl;
#endif

  permVec.resize(K_rows);
  SEQ_VECTOR < SEQ_VECTOR < eslocal > > 
    vec_I1_i2(K_rows, SEQ_VECTOR<eslocal >(2, 1));

  SEQ_VECTOR <eslocal > tmp_vec_s;
  tmp_vec_s.resize(sc_size);
  eslocal v1, n_mv, cnt_permut_vec;
  SEQ_VECTOR <eslocal >::iterator it;
  SEQ_VECTOR <eslocal > fix_dofs;
  fix_dofs.resize(sc_size);


  if (permutVectorActive<3){
    // set row permVec = {0,1,2,3,4,...,K_rows};
    if (fixing_nodes_or_dof==0 || (permutVectorActive==0)){
      for (eslocal i=0; i<K_rows; ++i) { permVec[i]=i;} // 0 1 2 K_rows-1
    }
    else
    {
      for (eslocal i=0; i<n_nodsSub; ++i) { permVec[i]=i;} // 0 1 2 n_nodsSub-1
    }
  }
//
  if (permutVectorActive==1){
    srand(time(NULL));
//    srand(0); // random will be constant until next compiling

    if (fixing_nodes_or_dof==0){
      random_shuffle ( permVec.begin(), permVec.end() );
    }
    else
    {
      std::srand(std::time(0));
      std::random_shuffle ( permVec.begin(), permVec.begin()+n_nodsSub);
      for (eslocal i=n_nodsSub;i>0;i--){
        for (eslocal j=0;j<dofPerNode;j++){
          permVec[dofPerNode*i-1-j] = dofPerNode*permVec[i-1]+j;
        }
      }
    }

    sort(permVec.begin(),permVec.begin()+nonsing_size);
    sort(permVec.begin()+nonsing_size,permVec.end());
  }
  else if (permutVectorActive==2){
    // random permutation
    n_mv = 0;                     // n_mv = size(unique(tmp_vec_s)) has to be equal to sc_size
    cnt_permut_vec=0;
    srand(time(NULL));
    // loop controls, if series 'tmp_vec_s' with unique integers has suffisciant dimension.
    // If not, missing numbers are added and checked again.
    do {
      for (eslocal i = 0;i<(sc_size-n_mv);i++){
        v1 = rand() % K_rows;
        tmp_vec_s[n_mv+i]=v1;
      }
      it=tmp_vec_s.begin();
      std::sort(tmp_vec_s.begin(), tmp_vec_s.end());
      it = std::unique(tmp_vec_s.begin(), tmp_vec_s.end());
      n_mv = distance(tmp_vec_s.begin(),it);
      cnt_permut_vec++;
   } while (n_mv != sc_size && cnt_permut_vec < 100);
    //
    eslocal ik=0,cnt_i=0;
    for (eslocal i = 0;i<(int)permVec.size();i++){
      if (i==tmp_vec_s[ik]){
        permVec[ik+nonsing_size]=tmp_vec_s[ik];
        ik++;
      }
      else{
        permVec[cnt_i]=i;
        cnt_i++;
      }
    }
  }

#if DBG > 1
  std::cout << "permVec     : ";
  for (auto& ip : permVec)
    std::cout << ip << ' ';
  std::cout<< '\n';
#endif




  SpMat K_modif;

  Eigen::Map<Eigen::VectorXi> perm_vec(permVec.data(),permVec.size());
  K_modif = m_spmatK.twistedBy(perm_vec.asPermutation().inverse());

  std::vector<T> trScaleMat(K_rows,T(0,0,0));
  double d_i(0);
  for (int dof = 0; dof < K_rows; dof++)
  {
    d_i = 0;
    if (diagonalScaling) d_i = 1.0 / sqrt(K_modif.coeff(dof,dof));
    trScaleMat[dof] = T(dof,dof,d_i);
  }

  SpMat ScaleMatrix(K_rows,K_rows);
  ScaleMatrix.setFromTriplets(trScaleMat.begin(),trScaleMat.end());

  K_modif = ScaleMatrix * K_modif * ScaleMatrix;



  for (eslocal i = 0;i<sc_size;i++)
    fix_dofs[i]=permVec[nonsing_size + i];

  Eigen::MatrixXd S_ss;
  SpMat K_rr, K_ss, K_rs, K_sr;


  K_rr = K_modif.block(0,0,nonsing_size,nonsing_size);
  K_rs = K_modif.block(0,nonsing_size,nonsing_size,sc_size);
  K_sr = K_modif.block(nonsing_size,0,sc_size,nonsing_size);
  K_ss = K_modif.block(nonsing_size,nonsing_size,sc_size,sc_size);


  Eigen::PardisoLU<SpMat >  Krr_solver;
  Krr_solver.compute(K_rr);


  Eigen::MatrixXd K_rsDense = K_rs.toDense();

  S_ss = K_ss - K_sr * (Krr_solver.solve(K_rsDense));

// IDENTIFICATIONS OF ZERO EIGENVALUES
  Eigen::MatrixXd U_, V_;
  Eigen::VectorXd diagS;
  linalg::svd0(S_ss,U_, diagS ,V_);


  eslocal defect_K_in = 0;// R_s_cols;
  double ratio;
  eslocal itMax = twenty < S_ss.rows() ? twenty : S_ss.rows();
 // for (eslocal i = itMax-1; i > 0;i--){
  for (eslocal i = S_ss.rows() - itMax; i < S_ss.rows()-1;i++){
    ratio = fabs(diagS(i+1)/diagS(i));
    if (ratio < jump_in_eigenvalues_alerting_singularity){
      defect_K_in=S_ss.rows() - i - 1;
      break;
    }
  }
  if (get_n_first_and_n_last_eigenvals_from_dense_S!=0){
    int i1i = get_n_first_and_n_last_eigenvals_from_dense_S;
    if (i1i>S_ss.rows()){i1i=S_ss.rows();}
    std::cout<<"eigenvals of S d{1:" << i1i << "} and d{" <<
         S_ss.rows()-get_n_first_and_n_last_eigenvals_from_dense_S+2 << ":"<< S_ss.rows()<< "}\n";

    for (eslocal i = 0 ; i < S_ss.rows(); i++){
      if (i < get_n_first_and_n_last_eigenvals_from_dense_S ||
            i > S_ss.rows()-get_n_first_and_n_last_eigenvals_from_dense_S){
        std::cout<< i+1 <<":\t"<< diagS[i] << "\n";
      }
    }
  }
// --------------- CREATING KERNEL R_s FOR SINGULAR PART (SCHUR COMPLEMENT)
  std::cout << "defect_K_in: " << defect_K_in << std::endl;

  Eigen::MatrixXd R_s = V_.rightCols(defect_K_in);


// --------------- CREATING KERNEL R_r FOR NON-SINGULAR PART
  Eigen::MatrixXd R_r(K_rr.rows(),defect_K_in);
  if (K_rr.cols() != 0){
    Eigen::MatrixXd KrsRs = K_rs * R_s;
    R_r = Krr_solver.solve(KrsRs);
  }

  Eigen::MatrixXd kerK(K_rows,defect_K_in);

  int R_r_rows = R_r.rows();
  for (eslocal j = 0; j < kerK.cols(); j++)
  {
    for (eslocal i = 0; i < R_r_rows; i++){
      kerK(permVec[i],j) = R_r(i,j) * ScaleMatrix.coeff(i,i);
    }
    for (eslocal i = 0; i < R_s.rows(); i++){
      kerK(permVec[i+R_r_rows],j) = 
        (-1) * R_s(i,j) * ScaleMatrix.coeff(i+R_r_rows,i+R_r_rows);
    }
  }

#if DBG > 2
  auto KR = K * kerK;
  std::cout << "|| K * R || "  << KR.norm() << std::endl;
#endif


  Eigen::FullPivHouseholderQR<Eigen::MatrixXd>
    qr(kerK.rows(), kerK.cols());
  qr.setThreshold(1e-15);
  qr.compute(kerK);
  //TODO not sure if it's really optimized.
  // Is internally Q.cols() == defect_K_in ???
  Eigen::MatrixXd Q =
    qr.matrixQ() * Eigen::MatrixXd::Identity(kerK.rows(), kerK.cols());


#if DBG > 2
  auto col_rank = qr.rank();
  std::cout << "Q.dim = " <<  Q.rows() << ' ' << Q.cols() << '\n';
  std::cout << "rank = " << col_rank << '\n';
  auto KQ = K * Q;
  std::cout << "|| K * Q || "  << KQ.norm() << std::endl;

#endif
  return Q;

}


std::vector <int> StiffnessMatrix::GetNullPivots(const Eigen::MatrixXd& R)
{


  int n_col_ = R.cols();
  int n_row_ = R.rows();


  std::vector<double> N(n_col_ * n_row_);

  int cnt(0);
  for (int col = 0; col < n_col_; col++){
    for (int row = 0; row < n_row_; row++){
      N[cnt++]=R(row,col);
    }
  }




  auto ij = [&](int i_, int j_){ return i_ + j_ * R.rows();};

  std::cout << "ij = " << ij(1,1) << std::endl;

  std::vector <double>::iterator  it;
  int I,colInd,rowInd;
  double *tmpV = new double[n_row_];
  int *_nul_piv = new int [n_row_];
  double pivot;
  int tmp_int;
  for (int i = 0;i < n_row_;i++){
      _nul_piv[i]=i;
  }
  for (int j=0;j<n_col_;j++){
    it = max_element(N.begin(),N.end()-j*n_row_,
        [](double i, double j){ return fabs(i)<=fabs(j); });

    I = it - N.begin();
    colInd = I/n_row_;
    rowInd = I-colInd*n_row_;
    for (int k=0;k<n_col_-j;k++){
      tmpV[k] = N[ij(n_row_-1-j,k)];
      N[ij(n_row_-1-j,k)] = N[ij(rowInd,k)];
      N[ij(rowInd,k)]= tmpV[k];
    }
    tmp_int = _nul_piv[rowInd];
    _nul_piv[rowInd] = _nul_piv[n_row_-1-j];
    _nul_piv[n_row_-1-j] = tmp_int;
    if (n_col_ - 1 - j != colInd) {
        memcpy( tmpV, &(N[ij(0,n_col_-1-j)]) , sizeof( double ) * n_row_);
        memcpy( &(N[ij(0, n_col_- 1 - j)]), &(N[ij(0, colInd)]), sizeof( double ) * n_row_);
        memcpy( &(N[ij(0,colInd)]),tmpV , sizeof( double ) * n_row_);
    }
    pivot = N[ij(n_row_-1-j,n_col_-1-j)];
    for (int J=0;J<n_col_-j-1;J++){
      for (int I=0;I<n_row_-j;I++){
        N[ij(I,J)] -= N[ij(I,n_col_-1-j)]*N[ij(n_row_-1-j,J)]/pivot;
      }
    }
  }
  std::vector<int> null_pivots(n_col_);
  for (int i = 0;i<n_col_;i++){
    null_pivots[i] = _nul_piv[n_row_-1-i];
  }
  sort(null_pivots.begin(),null_pivots.end());

  delete [] _nul_piv;
  delete [] tmpV;

  return null_pivots;

}

void StiffnessMatrix::ApplyDirichletBC(SpMat& spmat, std::vector<int> dirInd){

  Eigen::VectorXd diagMat = spmat.diagonal();
  double meanVal = diagMat.sum() / diagMat.size();

  for (int iD = 0; iD < static_cast<int>(dirInd.size()); iD++)
  {
    int kD = dirInd[iD];
    for (SpMat::InnerIterator it(spmat, kD); it; ++it)
    {
      if (it.row() == it.col())
        it.valueRef() = meanVal;
      else
        it.valueRef() = 0;
    }
  }
}

void StiffnessMatrix::SetDirichletPrecond(
    std::vector<int>& I_DOFs, std::vector<int>& B_DOFs)
{


  if (m_preconditionerType!= DIRICHLET) return;

  m_B_DOFs = B_DOFs;

  SpMat Kcopy = m_spmatK;
  int K_rows = Kcopy.rows();

//  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P(K_rows);

  // out_mat = mat.twistedBy( perm_vec.asPermutation().inverse() );

  Eigen::VectorXi perm_vec(K_rows);

  int cnt(0);
  for (auto& ii : I_DOFs) perm_vec(cnt++) = ii;
  for (auto& ii : B_DOFs) perm_vec(cnt++) = ii;

  Kcopy = m_spmatK.twistedBy(perm_vec.asPermutation().inverse());

  int nI = (int) I_DOFs.size();
  int nB = (int) B_DOFs.size();
  


  SpMat K_II = Kcopy.block(0,0,nI,nI);
  m_K_IB = Kcopy.block(0,nI,nI,nB);
  m_K_BI = Kcopy.block(nI,0,nB,nI);
  m_K_BB = Kcopy.block(nI,nI,nB,nB);


  m_pardisoSolverDirichlet_KII.analyzePattern(K_II);
  m_pardisoSolverDirichlet_KII.factorize(K_II);

#if DBG > 4
  m_dbg_printStiffnessMatrix(K_II,"K_II");
  m_dbg_printStiffnessMatrix(m_K_IB,"K_IB");
  m_dbg_printStiffnessMatrix(m_K_BI,"K_BI");
  m_dbg_printStiffnessMatrix(m_K_BB,"K_BB");
  m_dbg_printStiffnessMatrix(m_spmatK,"m_spmatK");
  m_dbg_printStiffnessMatrix(Kcopy,"Kcopy");
#endif





}


// dbg##############################################
void StiffnessMatrix::m_dbg_printStiffnessMatrix(const SpMat& spmat,
    std::string name)
{
  if (name.size()==0) name = "matK";
  std::string fname = name + std::to_string(m_myrank) + ".txt";
  tools::printMatrix(spmat,fname);
}

void StiffnessMatrix::m_dbg_printStiffnessMatrixSingularValues()
{
  std::cout << "singular values of K: \n";
  auto singVals = linalg::svd0(m_spmatK);
  for (int iv = 0; iv < singVals.size() ; iv++)
    std::cout << singVals[iv] << ' ';
  std::cout << '\n';
}

