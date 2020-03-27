#include "../include/directSolver.hpp"


using namespace Eigen;
using namespace std;

 
void DirectSolver::sym_factor(SpMat &A)
{
 
	//maxfct, mnum, phase, error, msglvl;
	msglvl = 0; 
 


	if (msglvl<0)
		msglvl = 0;
	mtype = 11;               
	nrows = static_cast<int>(A.rows());
	nnz = static_cast<int>(A.nonZeros());
	MKL_INT i;
	nrhs = 1;               /* Number of right hand sides. */
	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		iparm[i] = 0;


	valuePtr.resize(nnz);
	for (int i = 0; i < nnz; i++){
		valuePtr[i] = A.valuePtr()[i]; 
 	}

	cout << endl;
	innerIndexPtr.resize(nnz);
	for (int i = 0; i < nnz; i++){
		innerIndexPtr[i] = A.innerIndexPtr()[i];
 	}
	cout << endl;
	outerIndexPtr.resize(nrows + 1);
	for (int i = 0; i < nrows + 1; i++){
		outerIndexPtr[i] = A.outerIndexPtr()[i];
 	}
 

	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[2] = 1;         /*   */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Not in use */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off
						  (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;        /* One- or zero-based indexing of columns and rows.*/
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	error = 0;            /* Initialize error flag */
	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++)
		pt[i] = 0;
	/* -------------------------------------------------------------------- */
	///* .. Reordering and Symbolic Factorization. This step also allocates */
	///* all memory that is necessary for the factorization. */
	///* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&nrows, 
		&valuePtr[0],
		&outerIndexPtr[0],
		&innerIndexPtr[0],
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	} 
}

void DirectSolver::num_factor()
{
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&nrows,
		&valuePtr[0],
		&outerIndexPtr[0],
		&innerIndexPtr[0],
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");

}

void DirectSolver::solve(VectorXd& b, VectorXd& x)
{

	phase = 33;
	iparm[7] = 2;
  x.resize(b.size());
	x.setZero();
//	vector<double> x_(b.size(), 0);
//	vector<double> b_(&b[0], b.data() + b.size()); 

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&nrows,
		&valuePtr[0],
		&outerIndexPtr[0],
		&innerIndexPtr[0],
		&idum, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);

	//from the std vector to the eigen vector
//	x = Map<VectorXd, Unaligned>(x_.data(), x_.size());


}

int DirectSolver::testUnsym2(){

	MatrixXd A(5, 5);
	A << 1., -1., 0., -3., 0.,
		-2., 5., 0., 0., 0.,
		0., 0., 4., 6., 4.,
		-4., 0., 2., 7., 0.,
		0., 8., 0., 0., -5.;

	
	SpMat Asp;
	Asp = A.sparseView();

	//Asp.makeCompressed();

	cout << "sparse(A): \n";
	cout << Asp << endl;

	sym_factor(Asp);
	num_factor();
	VectorXd b(5), x(5);
	b.setOnes();

	solve(b, x);

	for (int i = 0; i < 5; i++) {
		printf("\n x [%d] = % f", i, x(i));
	}
	printf("\n");
	cout << A;

	cout << endl;

	//VectorXd x2 = A.lu().solve(b);
	//for (int i = 0; i < 5; i++) {
	//	printf("\n x [%d] = % f", i, x2(i));
	//}
	//printf("\n");

	return 0;
}


int DirectSolver::testUnsym(){
	/*
	********************************************************************************
	*                              INTEL CONFIDENTIAL
	*   Copyright(C) 2004-2008 Intel Corporation. All Rights Reserved.
	*   The source code contained  or  described herein and all documents related to
	*   the source code ("Material") are owned by Intel Corporation or its suppliers
	*   or licensors.  Title to the  Material remains with  Intel Corporation or its
	*   suppliers and licensors. The Material contains trade secrets and proprietary
	*   and  confidential  information of  Intel or its suppliers and licensors. The
	*   Material  is  protected  by  worldwide  copyright  and trade secret laws and
	*   treaty  provisions. No part of the Material may be used, copied, reproduced,
	*   modified, published, uploaded, posted, transmitted, distributed or disclosed
	*   in any way without Intel's prior express written permission.
	*   No license  under any  patent, copyright, trade secret or other intellectual
	*   property right is granted to or conferred upon you by disclosure or delivery
	*   of the Materials,  either expressly, by implication, inducement, estoppel or
	*   otherwise.  Any  license  under  such  intellectual property  rights must be
	*   express and approved by Intel in writing.
	*
	********************************************************************************
	*   Content : MKL PARDISO C example
	*
	********************************************************************************
	*/
	/* -------------------------------------------------------------------- */
	/* Example program to show the use of the "PARDISO" routine */
	/* on symmetric linear systems */
	/* -------------------------------------------------------------------- */
	/* This program can be downloaded from the following site: */
	/* www.pardiso-project.org */
	/* */
	/* (C) Olaf Schenk, Department of Computer Science, */
	/* University of Basel, Switzerland. */
	/* Email: olaf.schenk@unibas.ch */
	/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
	/* PARDISO prototype. */

/*
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
*/
/*
#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif
*/ 
#define MKL_INT int
 
	/* Matrix data. */
	MKL_INT n = 5;
	MKL_INT ia[6] = { 1, 4, 6, 9, 12, 14 };
	MKL_INT ja[13] = { 1, 2, 4,
		1, 2,
		3, 4, 5,
		1, 3, 4,
		2, 5 };
	double a[13] = { 1.0, -1.0, -3.0,
		-2.0, 5.0,
		4.0, 6.0, 4.0,
		-4.0, 2.0, 7.0,
		8.0, -5.0 };
	MKL_INT mtype = 11; /* Real unsymmetric matrix */
	/* RHS and solution vectors. */
	double b[5], x[5];
	MKL_INT nrhs = 1; /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */
	/* -------------------------------------------------------------------- */
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
	/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 1; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
	/* -------------------------------------------------------------------- */
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
	/* Set right hand side to one. */
	for (i = 0; i < n; i++) {
		b[i] = 1;
	}
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a, ia, ja, &idum, &nrhs,
		iparm, &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");
	for (i = 0; i < n; i++) {
		printf("\n x [%d] = % f", i, x[i]);
	}
	printf("\n");
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	return 0;
}


DirectSolver::~DirectSolver(){

	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&nrows,
		&valuePtr[0],
		&outerIndexPtr[0],
		&innerIndexPtr[0],
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
}

