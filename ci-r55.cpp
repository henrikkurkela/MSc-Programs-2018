/*
Helium CI Program
Developed using Microsoft Visual C++ 2017

Required 3rd party libraries:

Boost 1.65.1 - Boost License 1.0 - Multiple Authors
Eigen 3.3.4 - MPL2 - Multiple Authors
valandil/wignerSymbols - LGPL 3.0 - Author Joey Dumont

Equations given in () refer to R. D. Cowan's book "The Theory of Atomic Structure and Spectra".
*/

#include <iostream>
#include <thread>
#include <mutex>
#include <cmath>
#include <string>
#include <chrono>
#include <boost/multi_array.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "wignerSymbols.h"

using namespace std;

/* --------------------------------------------------------------- */
/* The factorial function */
double factorial(int n) {
	if (n == 1 || n == 0) return 1;
	else return n * factorial(n - 1);
}
/* The Kronecker Delta function */
int kdelta(int a, int b) {
	if (a == b) return 1;
	else return 0;
}
/* The ck function (11.23) */
double ck(int l, int ld, int k) {
	double rval = 0;
	rval = pow(-1.00, double(l)) * pow(double(((2 * l + 1) * (2 * ld + 1))), 0.5) * WignerSymbols::wigner3j(l, k, ld, 0, 0, 0);
	return rval;
}
/* The alpha function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double alfa(int n, int l, int x, double z) {
	double rval = 0;
	/* Calculate the alpha(n,l,x) and place the result in rval */
	rval = 
		pow(-1, x) *
		pow(2, l + x + 1) *
		pow(z, l + x + 1.5) *
		pow(pow(n, -1), l + x + 2) *
		pow(factorial(n - l - 1) * factorial(n + l), 0.5)
		/ (factorial(n - l - 1 - x) * factorial(2 * l + 1 + x) * factorial(x));
	return rval;
}
/* The single-electron Hamiltonian matrix element function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double hamiltonian(int ni, int li, int mi, double zi, int nj, int lj, int mj, double zj, int z) {
	double omegaij = (ni * nj) / (ni * zj + nj * zi);
	double rval = 0;
	if (li != lj) return 0;
	if (mi != mj) return 0;
	for (int xi = 0; xi < ni - li; xi++) {
		for (int xj = 0; xj < nj - lj; xj++) {
			rval +=
				alfa(ni, li, xi, zi) * alfa(nj, lj, xj, zj) * pow(omegaij, 2 * li + xi + xj + 1) * (
					z * omegaij * factorial(2 * li + xi + xj + 1) + 0.5 * factorial(2 * li + xi + xj) * (
					(zj * omegaij / nj) * (2 * li + xi + xj + 1) * (
						(zj * omegaij / nj) * (2 * li + xi + xj + 2) - 2 * (li + xj + 1)
						) + xj * (2 * li + xj + 1)
						)
					);
		}
	}
	return -rval;
}
/* The rkijtu function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double rkijtu(int ni, int li, double zi, int nj, int lj, double zj, int nt, int lt, double zt, int nu, int lu, double zu, int k) {
	double omegait = (ni * nt) / (ni * zt + nt * zi);
	double omegaju = (nj * nu) / (nj * zu + nu * zj);
	double rval = 0;
	double temp = 0;
	for (int xi = 0; xi < ni - li; xi++) {
		for (int xj = 0; xj < nj - lj; xj++) {
			for (int xt = 0; xt < nt - lt; xt++) {
				for (int xu = 0; xu < nu - lu; xu++) {
					temp = 0;
					for (int b = 0; b < lj + lu + xj + xu + k + 3; b++) {
						temp += 1.00 / (factorial((b)) * pow(omegaju, b)) * pow((1.00 / omegaju) + (1.00 / omegait), k - li - lt - xi - xt - b - 2) * factorial((b + li + lt + xi + xt - k + 1));
					}
					rval += alfa(ni, li, xi, zi) * alfa(nj, lj, xj, zj) * alfa(nt, lt, xt, zt) * alfa(nu, lu, xu, zu) * pow(omegaju, k + lj + lu + xj + xu + 3) * factorial((lj + lu + xj + xu + k + 2)) *
						(pow(omegait, li + lt + xi + xt - k + 2) * factorial(li + lt + xi + xt - k + 1) - temp);
				}
			}
		}
	}
	return double(rval);
}
/* The Rkijtu function combining two rkijtu's to create (6.23) */
double Rkijtu(int ni, int li, int nj, int lj, int nt, int lt, int nu, int lu, double z, int k) {
	double rval = 0;
	rval = rkijtu(ni, li, z, nj, lj, z, nt, lt, z, nu, lu, z, k) + rkijtu(nj, lj, z, ni, li, z, nu, lu, z, nt, lt, z, k);
	return rval;
}
/* --------------------------------------------------------------- */
/* The direct coefficient (13.22) */
double rkd(int l, int s, int ld, int sd, int li, int lj, int lt, int lu, int k) {
	double rval = 0;
	rval =
		kdelta(l, ld) * kdelta(s, sd) * pow(-1, double(lj + lt + l)) *
		ck(li, lt, k) * ck(lj, lu, k) * WignerSymbols::wigner6j(li, lj, l, lu, lt, k);
	return rval;
}
/* The exchange coefficient (13.23) */
double rke(int l, int s, int ld, int sd, int li, int lj, int lt, int lu, int k) {
	double rval = 0;
	rval =
		kdelta(l, ld) * kdelta(s, sd) * pow(-1, double(s)) *
		ck(li, lu, k) * ck(lt, lj, k) * WignerSymbols::wigner6j(li, lj, l, lt, lu, k);
	return rval;
}
/* The Coulomb matrix element function (13.21) */
double coulomb(int l, int s, int ld, int sd, int ni, int li, int nj, int lj, int nt, int lt, int nu, int lu, double z) {
	double rval = 0;
	int kmin = -INT_MAX;
	int kmax = INT_MAX;
	/* Figure out the k values for the direct part */
	if (kmin < abs(li - lt)) kmin = abs(li - lt);
	if (kmin < abs(lj - lu)) kmin = abs(lj - lu);
	if (kmax > (li + lt)) kmax = (li + lt);
	if (kmax > (lj + lu)) kmax = (lj + lu);
	if ((kmax % 2) == (kmin % 2))
		for (int k = kmin; k <= kmax; k += 2) {
			rval = rval +
				(0.5 * rkd(l, s, ld, sd, li, lj, lu, lt, k) * Rkijtu(ni, li, nj, lj, nt, lt, nu, lu, z, k));
		};
	kmin = -INT_MAX;
	kmax = INT_MAX;
	/* Figure out the k values for the exchange part */
	if (kmin < abs(li - lu)) kmin = abs(li - lu);
	if (kmin < abs(lj - lt)) kmin = abs(lj - lt);
	if (kmax > (li + lu)) kmax = (li + lu);
	if (kmax > (lj + lt)) kmax = (lj + lt);

	if ((kmax % 2) == (kmin % 2))
		for (int k = kmin; k <= kmax; k += 2) {
			rval = rval +
				(0.5 * rke(l, s, ld, sd, li, lj, lu, lt, k) * Rkijtu(ni, li, nj, lj, nu, lu, nt, lt, z, k));
		};
	return rval;
}
/* --------------------------------------------------------------- */

/* A structure defining a two-electron state. */
struct csf_twoelectron {
	int l;
	int s;
	int n1;
	int l1;
	int n2;
	int l2;
};

/* The Helium basis creation function. Looks up pairs of electrons that have the required LS term. */
int createbasis_helium(int l, int s, int nmax, int lmax, csf_twoelectron **basis) {
	/* rval will be the total size of the basis, and will be returned */
	int rval = 0;
	int index = 0;
	int parity = 0;
	csf_twoelectron *pointer = NULL;
	for (int n = 1; n < nmax + 1; n++) {
		for (int nd = n; nd < nmax + 1; nd++) {
			for (int ml = 0; ml < min(lmax + 1, n); ml++) {
				for (int mld = ml; mld < min(lmax + 1, nd); mld++) {
					/* Parity check */
					if ((ml + mld) % 2 == parity) {
						/* Now we have two cases, either the electrons are equivalent, or not.
						If electrons are non-equivalent, we can skip the S check because for two
						electrons there always exists a singlet and a triplet because individual spins are not restricted.
						If the electrons are equivalent, only singlet state can exist.
						*/
						/* The case of singlet configuration. No equivalency check done. */
						if (s == 0) {
							for (int i = abs(ml - mld); i < ml + mld + 1; i++) {
								if (i == l) rval++;
							}
						}
						/* The case of triplet configuration. Check for equivalent electrons and reject them. */
						else if (s == 1) {
							for (int i = abs(ml - mld); i < ml + mld + 1; i++) {
								if ((i == l) && ((n != nd) || (ml != mld))) rval++;
							}
						}
						/* Two electrons can only have singlets and triplets. If S is neither 0 or 1, reject all */
						else {
							cout << "Spin violation. S must be either 0 or 1 for two electron configurations." << endl;
							return -1;
						}
					}
				}
			}
		}
	}

	cout << "Total configuration basis size: " << rval << endl;
	/* Reassign the given csf struct pointer, and copy the address for later use */
	*basis = new csf_twoelectron[rval];
	pointer = *basis;
	for (int n = 1; n < nmax + 1; n++) {
		for (int nd = n; nd < nmax + 1; nd++) {
			for (int ml = 0; ml < min(lmax + 1, n); ml++) {
				for (int mld = ml; mld < min(lmax + 1, nd); mld++) {
					/* Parity check */
					if ((ml + mld) % 2 == parity) {
						/* The case of singlet configuration. No equivalency check done. */
						if (s == 0) {
							for (int i = abs(ml - mld); i < ml + mld + 1; i++) {
								if (i == l) {
									pointer[index].l = l;
									pointer[index].s = s;
									pointer[index].n1 = n;
									pointer[index].l1 = ml;
									pointer[index].n2 = nd;
									pointer[index].l2 = mld;
									index++;
								}
							}
						}
						/* The case of triplet configuration. Check for equivalent electrons and reject them. */
						else if (s == 1) {
							for (int i = abs(ml - mld); i < ml + mld + 1; i++) {
								if ((i == l) && ((n != nd) || (ml != mld))) {
									pointer[index].l = l;
									pointer[index].s = s;
									pointer[index].n1 = n;
									pointer[index].l1 = ml;
									pointer[index].n2 = nd;
									pointer[index].l2 = mld;
									index++;
								}
							}
						}
					}
				}
			}
		}
	}

	return rval;
}

/*
Default values for various inputs. Some will be altered on execution.
LSL is the LS term L (this can not be overridden)
LSS is the LS term S (this can not be overridden)
Z is the atomic number (this can not be overridden)
ZEFF is the effective nuclear charge
MAXN is the maximum N quantum number used to create the configuration basis
MAXL is the maximum L quantum number used to create the configuration basis
*/
#define LSL 0
#define LSS 0
#define Z 2
#define ZEFF 2.000
#define MAXN 3
#define MAXL 2

/* The calculation variables. Some of these may be altered by user input. */
int csfsize = 1;
int l = LSL;
int s = LSS;
int maxn = MAXN;
int maxl = MAXL;
int workers = thread::hardware_concurrency();
double zeff = ZEFF;
int workqueue[2] = { 0, 0 };
mutex workqueue_mutex;
int reserved_elements = 0;
int total_elements = 0;
/* Pointer variables. These should be deleted before exit. */
csf_twoelectron *csfs;
thread *threads;
/* Eigen matrices and solvers */
Eigen::MatrixXd hmtx;
Eigen::MatrixXd cmtx;
Eigen::MatrixXd fmtx;
Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
/* The threadrunner function that calculates the Coulomb matrix elements. */
int threadrunner(int index) {
	int i;
	int j;
	/* Figure out the lowest non-calculated coulomb matrix item */
	while (true) {
		{
			lock_guard<mutex> guard(workqueue_mutex);
			/* If within bounds, reserve the matrix element */
			if ((workqueue[0] < csfsize) && (workqueue[1] < csfsize)) {
				i = workqueue[0];
				j = workqueue[1];
				reserved_elements += 2;
			}
			else {
				return 0;
			}
			/* Increase the column value */
			workqueue[1]++;
			/* Check for row advance */
			if (workqueue[1] > csfsize - 1) {
				workqueue[0]++;
				workqueue[1] = workqueue[0];
			}
		}
		/* Calculate the reserved element */
		cmtx(i, j) = coulomb(csfs[i].l, csfs[i].s, csfs[j].l, csfs[j].s, csfs[i].n1, csfs[i].l1, csfs[i].n2, csfs[i].l2, csfs[j].n1, csfs[j].l1, csfs[j].n2, csfs[j].l2, zeff);
		cmtx(j, i) = cmtx(i, j);
	}
	return 0;
}

/* The main program loop */
int main(int argc, char * argv[]) {
	cout << "******** CI calculation setup *********" << endl
		 << "Give the amount of threads to run: ";
	cin >> workers;
	cout << "Give the effective nuclear charge Z: ";
	cin >> zeff;
	cout << "Give the maximum value for N: ";
	cin >> maxn;
	cout << "Give the maximum value for L: ";
	cin >> maxl;
	/* Get the current time */
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	/* Set the thread vector length to match specified thread count */
	threads = new thread[workers];
	/* Call the createbasis -function with the LS term and maximum values for N and L */
	csfsize = createbasis_helium(l, s, maxn, maxl, &csfs);
	/* Resize the Hamiltonian and Coulomb matrices to configuration state function basis size */
	hmtx.resize(csfsize, csfsize);
	cmtx.resize(csfsize, csfsize);
	fmtx.resize(csfsize, csfsize);
	total_elements = csfsize * csfsize;
	cout << "********** Energy level list **********" << endl;
	/* Loop through all configurations and calculate single-electron Hamiltonians. */
	for (int i = 0; i < csfsize; i++) {
		for (int j = i; j < csfsize; j++) {
			/* Check for inequivalent electrons */
			bool lhs_equal = false;
			bool rhs_equal = false;
			if ((csfs[i].n1 == csfs[i].n2) && (csfs[i].l1 == csfs[i].l2)) lhs_equal = true;
			if ((csfs[j].n1 == csfs[j].n2) && (csfs[j].l1 == csfs[j].l2)) rhs_equal = true;
			/* Select the correct case as given in "Overlap integrals in LS-coupling" by K. Jänkälä */
			if ((rhs_equal == true) && (lhs_equal == true)) {
				hmtx(i, j) = 2 * hamiltonian(csfs[i].n1, csfs[i].l1, 0, zeff, csfs[j].n1, csfs[j].l1, 0, zeff, Z) * kdelta(csfs[i].n1, csfs[j].n1) * kdelta(csfs[i].l1, csfs[j].l1);
			}
			else if ((rhs_equal == false) && (lhs_equal == false)) {
				hmtx(i, j) =
					hamiltonian(csfs[i].n1, csfs[i].l1, 0, zeff, csfs[j].n1, csfs[j].l1, 0, zeff, Z) * kdelta(csfs[i].n2, csfs[j].n2) * kdelta(csfs[i].l2, csfs[j].l2) +
					hamiltonian(csfs[i].n2, csfs[i].l2, 0, zeff, csfs[j].n2, csfs[j].l2, 0, zeff, Z) * kdelta(csfs[i].n1, csfs[j].n1) * kdelta(csfs[i].l1, csfs[j].l1) +
					hamiltonian(csfs[i].n1, csfs[i].l1, 0, zeff, csfs[j].n2, csfs[j].l2, 0, zeff, Z) * kdelta(csfs[i].n1, csfs[j].n1) * kdelta(csfs[i].l1, csfs[j].l1) +
					hamiltonian(csfs[i].n2, csfs[i].l2, 0, zeff, csfs[j].n1, csfs[j].l1, 0, zeff, Z) * kdelta(csfs[i].n1, csfs[j].n2) * kdelta(csfs[i].l1, csfs[j].l2);
			}
			else {
				hmtx(i, j) =
					pow(2, 0.5) * hamiltonian(csfs[i].n1, csfs[i].l1, 0, zeff, csfs[j].n1, csfs[j].l1, 0, zeff, Z) * kdelta(csfs[i].n2, csfs[j].n2) * kdelta(csfs[i].l2, csfs[j].l2) +
					pow(2, 0.5) * hamiltonian(csfs[i].n2, csfs[i].l2, 0, zeff, csfs[j].n2, csfs[j].l2, 0, zeff, Z) * kdelta(csfs[i].n1, csfs[j].n1) * kdelta(csfs[i].l1, csfs[j].l1);
			}
			if (i != j) hmtx(j, i) = hmtx(i, j);
		}
	}
	/* Spawn threads to calculate the Coulomb matrix elements. */
	for (int i = 0; i < workers; i++) {
		threads[i] = thread(threadrunner, i);
	}
	/* Join all threads */
	for (int i = 0; i < workers; i++) {
		threads[i].join();
	}
	fmtx = hmtx + cmtx;
	/* Use the eigensolver to solve the eigenvalue problem (13.9) */
	solver.compute(fmtx, Eigen::MatrixXd::Identity(csfsize, csfsize));
	cout << solver.eigenvalues() << endl;
	/* Display the total calculation time */
	std::chrono::system_clock::time_point stop = std::chrono::system_clock::now();
	std::chrono::duration<double> runtime = stop - start;
	cout << "Runtime: " << runtime.count() << " s.";
	/* Delete allocated pointer variables and exit */
	delete[] csfs;
	delete[] threads;
	return 0;
}