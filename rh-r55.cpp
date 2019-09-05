/*
Roothaan-Hall Program Using Symbolic Matrix Element Integrals
Developed using Microsoft Visual Studio 2017

Required 3rd party libraries:

Boost 1.65.1 - Boost License 1.0 - Multiple Authors
Eigen 3.3.4 - MPL2 - Multiple Authors
valandil/wignerSymbols - LGPL 3.0 - Author Joey Dumont

Equations in parentheses refer to R. D. Cowan's book "The Theory of Atomic Structure and Spectra".

Acknowledgements:
https://www.eriksmistad.no/measuring-runtime-in-milliseconds-using-the-c-11-chrono-library/ - Usage of the standard chrono library
*/

#include <iostream>
#include <thread>
#include <chrono>
#include <mutex>
#include <cmath>
#include <string>
#include <boost/multi_array.hpp>
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
/* The Ck function (11.23) */
double ckq(int l, int ld, int m, int md, int k, int q) {
	double rval = 0;
	rval = pow(-1.00, -m) * pow(double(((2 * l + 1) * (2 * ld + 1))), 0.5) * WignerSymbols::wigner3j(l, k, ld, 0, 0, 0) * WignerSymbols::wigner3j(l, k, ld, -m, q, md);
	return rval;
}
/* The alpha function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double alfa(int n, int l, int x, double z) {
	/* Declare a Boost 100-decimal floating point variable rval */
	double rval = 0;
	/* Calculate the alpha(n,l,x) and place the result in rval */
	rval = (
		pow(-1, x) *
		pow(2, l + x + 1) *
		pow(z, l + x + 1.5) *
		pow(pow(n, -1), l + x + 2) *
		pow(double(factorial(n - l - 1)) * double(factorial(n + l)), 0.5)
		/ (double(factorial(n - l - 1 - x)) * double(factorial(2 * l + 1 + x)) * double(factorial(x)))
		);
	return rval;
}
/* The overlap matrix element function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double overlap(int ni, int li, int mi, double zi, int nj, int lj, int mj, double zj) {
	double omegaij = (ni * nj) / (ni * zj + nj * zi);
	double rval = 0;
	/* Check the L and M quantum numbers. If they do not match, the orbitals are not overlapping. */
	if (li != lj) return 0;
	if (mi != mj) return 0;
	for (int xi = 0; xi < ni - li; xi++) {
		for (int xj = 0; xj < nj - lj; xj++) {
			rval += alfa(ni, li, xi, zi) * alfa(nj, lj, xj, zj) * pow(omegaij, 2 * li + xi + xj + 3) * factorial(2 * li + xi + xj + 2);
		}
	}
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
	return - rval;
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
						temp += 1 / (factorial(b) * pow(omegaju, b)) * pow((1.00 / omegaju) + (1.00 / omegait), k - li - lt - xi - xt - b - 2) * factorial(b + li + lt + xi + xt - k + 1);
					}
					rval += alfa(ni, li, xi, zi) * alfa(nj, lj, xj, zj) * alfa(nt, lt, xt, zt) * alfa(nu, lu, xu, zu) * pow(omegaju, k + lj + lu + xj + xu + 3) * factorial(lj + lu + xj + xu + k + 2) *
						(pow(omegait, li + lt + xi + xt - k + 2) * factorial(li + lt + xi + xt - k + 1) - temp);
				}
			}
		}
	}
	return rval;
}
/* The Coulomb matrix element function given in "Analytic matrix elements for screened hydrogenic functions" by K. Jänkälä */
double hijtu(int ni, int li, int mi, double zi, int nj, int lj, int mj, double zj, int nt, int lt, int mt, double zt, int nu, int lu, int mu, double zu) {
	double rval = 0;
	int mink = max(abs(li - lt), abs(lj - lu));
	int maxk = min(li + lt, lj + lu);

	for (int k = mink; k < maxk + 1; k += 2) {
		for (int q = -k; q < k + 1; q++) {
			rval += pow(-1, q) * (rkijtu(ni, li, zi, nj, lj, zj, nt, lt, zt, nu, lu, zu, k) + rkijtu(nj, lj, zj, ni, li, zi, nu, lu, zu, nt, lt, zt, k)) *
				ckq(li, lt, mi, mt, k, -q) * ckq(lj, lu, mj, mu, k, q);
		}
	}
	return rval;
}
/* --------------------------------------------------------------- */

/* 
A structure used in solving Roothaan-Hall calculation. 
Holds information about basis functions and calculated solution.
*/
struct dataset {
	int size = 1;
	int z = 2;
	boost::multi_array<int, 1> nbase;
	boost::multi_array<int, 1> lbase;
	boost::multi_array<int, 1> mbase;
	boost::multi_array<double, 1> zbase;
	boost::multi_array<double, 4> coulomb;
	Eigen::VectorXd lowest_eigenvector;
	Eigen::MatrixXd all_eigenvectors;
	double lowest_eigenvalue = 0;
	double total_electron_energy = 0;
};

/* Main Roothaan-Hall calculation routine */
double calculation_runner(dataset &data) {
	double rval = 0;
	double rval_old = 0;
	double lowest_eigenvalue = INT_MAX;
	int lowest_orbital = 0;
	int iterations = 100;
	int orbital_occupation_count = data.z / 2;
	int *occupation_order = new int[orbital_occupation_count];
	double min_convergence = 0.000000001;

	Eigen::MatrixXd overlap_matrix(data.size, data.size);
	Eigen::MatrixXd hamiltonian_matrix(data.size, data.size);
	Eigen::MatrixXd coefficient_matrix(data.size, data.size);
	Eigen::MatrixXd density_matrix(data.size, data.size);
	Eigen::MatrixXd fock_matrix(data.size, data.size);
	Eigen::VectorXd norm_vector(data.size);
	Eigen::VectorXd eigen_values(data.size);
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;

	for (int i = 0; i < data.size; i++) {
		for (int j = 0; j < data.size; j++) {
			overlap_matrix(i, j) = overlap(data.nbase[i], data.lbase[i], data.mbase[i], data.zbase[i], data.nbase[j], data.lbase[j], data.mbase[j], data.zbase[j]);
			hamiltonian_matrix(i, j) = hamiltonian(data.nbase[i], data.lbase[i], data.mbase[i], data.zbase[i], data.nbase[j], data.lbase[j], data.mbase[j], data.zbase[j], data.z);
			for (int t = 0; t < data.size; t++) {
				for (int u = 0; u < data.size; u++) {
					data.coulomb[i][j][t][u] = hijtu(
						data.nbase[i], data.lbase[i], data.mbase[i], data.zbase[i], 
						data.nbase[j], data.lbase[j], data.mbase[j], data.zbase[j], 
						data.nbase[t], data.lbase[t], data.mbase[t], data.zbase[t], 
						data.nbase[u], data.lbase[u], data.mbase[u], data.zbase[u]);
				}
			}
		}
	}
	/* Set up the iteration loop by solving with single-electron Hamiltonian */
	solver.compute(hamiltonian_matrix, overlap_matrix);
	coefficient_matrix = solver.eigenvectors().real();
	/* Find the initial lowest orbital ordering scheme */
	eigen_values = solver.eigenvalues().real();
	for (int i = 0; i < orbital_occupation_count; i++) {
		lowest_eigenvalue = INT_MAX;
		for (int j = 0; j < data.size; j++) {
			if (lowest_eigenvalue > eigen_values(j)) {
				lowest_eigenvalue = eigen_values(j);
				occupation_order[i] = j;
			}
		}
		eigen_values(occupation_order[i]) = INT_MAX;
	}
	/* Renormalize the initial eigenvectors */
	for (int k = 0; k < data.size; k++) {
		norm_vector(k) = 0;
		for (int i = 0; i < data.size; i++) {
			norm_vector(k) += pow(coefficient_matrix(i, k), 2);
		}
		for (int j = 0; j < data.size - 1; j++) {
			for (int jd = j + 1; jd < data.size; jd++) {
				norm_vector(k) += 2 *
					coefficient_matrix(j, k) *
					coefficient_matrix(jd, k) *
					overlap_matrix(j, jd);
			}
		}
	}
	for (int i = 0; i < data.size; i++) {
		for (int j = 0; j < data.size; j++) {
			coefficient_matrix(i, j) = pow(norm_vector(j), -0.5) * coefficient_matrix(i, j);
		}
	}
	/* Create the initial density matrix */
	for (int i = 0; i < data.size; i++) {
		for (int j = 0; j < data.size; j++) {
			density_matrix(i, j) = 0;
			for (int k = 0; k < orbital_occupation_count; k++) {
				density_matrix(i, j) += 2 * coefficient_matrix(i, occupation_order[k]) * coefficient_matrix(j, occupation_order[k]);
			}
		}
	}
	/* Self-Consistent Field iteration loop */
	for (int n = 0; n < iterations; n++) {
		/* Reevaluate the Fock Hamiltonian */
		for (int i = 0; i < data.size; i++) {
			for (int t = 0; t < data.size; t++) {
				fock_matrix(i, t) = hamiltonian_matrix(i, t);
				for (int j = 0; j < data.size; j++) {
					for (int u = 0; u < data.size; u++) {
						fock_matrix(i, t) += density_matrix(j, u) * 
							(data.coulomb[i][j][t][u] - 0.5 * data.coulomb[i][j][u][t]);
					}
				}
			}
		}
		/* Solve the eigenvalue problem */
		solver.compute(fock_matrix, overlap_matrix);
		coefficient_matrix = solver.eigenvectors().real();
		/* Find the lowest orbital */

		/* Find the lowest orbital ordering scheme */
		eigen_values = solver.eigenvalues().real();
		for (int i = 0; i < orbital_occupation_count; i++) {
			lowest_eigenvalue = INT_MAX;
			for (int j = 0; j < data.size; j++) {
				if (lowest_eigenvalue > eigen_values(j)) {
					lowest_eigenvalue = eigen_values(j);
					occupation_order[i] = j;
				}
			}
			eigen_values(occupation_order[i]) = INT_MAX;
		}
		/* Renormalize the eigenvectors */
		for (int k = 0; k < data.size; k++) {
			norm_vector(k) = 0;
			for (int i = 0; i < data.size; i++) {
				norm_vector(k) += pow(coefficient_matrix(i, k), 2.00);
			}
			for (int j = 0; j < data.size - 1; j++) {
				for (int jd = j + 1; jd < data.size; jd++) {
					norm_vector(k) += 2 *
						coefficient_matrix(j, k) *
						coefficient_matrix(jd, k) *
						overlap_matrix(j, jd);
				}
			}
		}
		for (int i = 0; i < data.size; i++) {
			for (int j = 0; j < data.size; j++) {
				coefficient_matrix(i, j) = pow(norm_vector(j), -0.5) * coefficient_matrix(i, j);
			}
		}
		/* Create initial density matrix */
		for (int i = 0; i < data.size; i++) {
			for (int j = 0; j < data.size; j++) {
				density_matrix(i, j) = 0;
				for (int k = 0; k < orbital_occupation_count; k++) {
					density_matrix(i, j) += 2 * coefficient_matrix(i, occupation_order[k]) * coefficient_matrix(j, occupation_order[k]);
				}
			}
		}
		/* Calculate the total electronic charge density */
		/* When the program functions well enough, this step can be removed but it provides good debugging insight (when printed). */
		double temp = 0;
		for (int a = 0; a < orbital_occupation_count; a++) {
			for (int i = 0; i < data.size; i++) {
				for (int j = 0; j < data.size; j++) {
					temp += 2 *
						coefficient_matrix(i, occupation_order[a]) *
						coefficient_matrix(j, occupation_order[a]) *
						overlap_matrix(i, j);
				}
			}
		}
		rval_old = rval;
		rval = 0;
		/* Calculate the total electronic energy */
		for (int i = 0; i < data.size; i++) {
			for (int j = 0; j < data.size; j++) {
				rval += density_matrix(i, j) * hamiltonian_matrix(i, j);
				for (int t = 0; t < data.size; t++) {
					for (int u = 0; u < data.size; u++) {
						rval += 0.5 *
							density_matrix(i, t) *
							density_matrix(j, u) *
							(data.coulomb[i][j][t][u] - 0.5 * data.coulomb[i][j][u][t]);
					}
				}
			}
		}
		/* Check convergence criterion */
		if (abs(rval_old - rval) < min_convergence) {
			n = iterations;
		}
		else if (n == iterations - 1) {
			data.lowest_eigenvector.resize(1);
			data.lowest_eigenvector(0) = 10;
			return 0;
		}
	}
	/* Copy lowest eigenvector and eigenvalue in to the data set */
	data.lowest_eigenvector.resize(data.size);
	data.lowest_eigenvector = coefficient_matrix.col(occupation_order[0]);
	for (int i = 0; i < data.size; i++) {
		if (solver.eigenvalues().real()(i) < data.lowest_eigenvalue) {
			data.lowest_eigenvalue = solver.eigenvalues().real()(i);
		}
	}
	/* Copy all linear expansion coefficients in to the data set */
	data.all_eigenvectors.resize(data.size, data.size);
	data.all_eigenvectors = coefficient_matrix;
	/* Copy the total electronic energy in to the data set */
	data.total_electron_energy = rval;
	/* Return the total electronic energy */
	return rval;
}

/* Global multi-threading variables */
mutex worker_mutex;
mutex cout_mutex;
int next_calculation = 0;
int last_calculation = 0;
int concurrency_number = 2;
dataset *table_element_data;
boost::multi_array<int, 1> nbase;
boost::multi_array<int, 1> lbase;
boost::multi_array<int, 1> mbase;
boost::multi_array<double, 1> zbase;

/* This function is used to drive Roothaan-Hall calculations when using table mode */
void thread_runner() {
	int index = 0;
	double rval = 0;
	while (true) {
		{
			lock_guard<mutex> lock(worker_mutex);
			if (next_calculation < last_calculation) {
				index = next_calculation;
				next_calculation++;
			}
			else return;
		}

		rval = calculation_runner(table_element_data[index]);
	}
	return;
}

/* This function is used to increment effective nuclear charge numbers between grid points when using table mode */
void increment_zbase(int index, int maxstep, double min, double increment) {
	if ( abs(zbase[index] - (min + maxstep * increment)) < 0.0001) {
		zbase[index] = min;
		if (index > 0) increment_zbase(index - 1, maxstep, min, increment);
	}
	else {
		zbase[index] += increment;
	}
	return;
}

/* This function is used to calculate the gradient at current point when using gradient mode */
void calculate_gradient(dataset &data, boost::multi_array<double, 1> &direction) {
	/* Create a local copy of the data set in argument */
	dataset gradient_data;
	gradient_data.z = data.z;
	gradient_data.size = data.size;
	gradient_data.nbase.resize(boost::extents[data.size]);
	gradient_data.nbase = data.nbase;
	gradient_data.lbase.resize(boost::extents[data.size]);
	gradient_data.lbase = data.lbase;
	gradient_data.mbase.resize(boost::extents[data.size]);
	gradient_data.mbase = data.mbase;
	gradient_data.zbase.resize(boost::extents[data.size]);
	gradient_data.zbase = data.zbase;
	gradient_data.coulomb.resize(boost::extents[data.size][data.size][data.size][data.size]);
	/* Calculate the gradient in all directions  */
	for (int i = 0; i < data.size; i++) {
		direction[i] = 0;
		direction[i] += calculation_runner(gradient_data);
		gradient_data.zbase[i] += 0.01;
		direction[i] -= calculation_runner(gradient_data);
		gradient_data.zbase[i] -= 0.01;
	}
	return;
}

/* The main program loop */
int main(int argc, char *argv[]) {
	/* Argument vector variables */
	bool table_calculation_mode = false;
	bool gradient_calculation_mode = false;
	bool enable_full_concurrency = false;
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-table") == 0) table_calculation_mode = true;
		if (strcmp(argv[i], "-gradient") == 0) gradient_calculation_mode = true;
		if (strcmp(argv[i], "-enable") == 0) enable_full_concurrency = true;
	}
	if (enable_full_concurrency == true) {
		cout << "[SETUP] Worker count set to 1 thread per 1 logical CPU!" << endl;
		concurrency_number = 1;
	}
	else {
		cout << "[SETUP] Assuming Hyper-Threaded CPU, worker count set to 1 thread per 2 logical CPU's!" << endl;
	}
	/* User has selected the table calculation mode */
	if (table_calculation_mode == true) {
		thread *worker;
		int functions;
		double min_z;
		double step_z;
		int steps;
		int z;
		int lowest_energy_index = 0;
		double lowest_energy = INT_MAX;
		double lowest_energy_sanity_check = -INT_MAX;
		auto start_time = chrono::high_resolution_clock::now();
		cout << "[SETUP] Table mode enabled!" << endl;
		cout << "Give the number of basis functions: ";
		cin >> functions;
		cout << "Give the element number (MUST BE EVEN): ";
		cin >> z;
		cout << "Give the minimum effective nuclear charge: ";
		cin >> min_z;
		cout << "Give the effective nuclear charge step size: ";
		cin >> step_z;
		cout << "Give the effective nuclear charge step count: ";
		cin >> steps;
		worker = new thread[thread::hardware_concurrency() / concurrency_number];
		nbase.resize(boost::extents[functions]);
		lbase.resize(boost::extents[functions]);
		mbase.resize(boost::extents[functions]);
		zbase.resize(boost::extents[functions]);
		for (int i = 0; i < functions; i++) {
			cout << "Give the N L M quantum numbers of basis function " << i + 1 << ": ";
			cin >> nbase[i];
			cin >> lbase[i];
			cin >> mbase[i];
			zbase[i] = min_z;
		}
		start_time = std::chrono::high_resolution_clock::now();
		table_element_data = new dataset[pow(steps + 1, functions)];
		cout << "Total calculation space size: " << pow(steps + 1, functions) << endl;
		cout << "Setting up calculation... " << endl;
		for (int i = 0; i < pow(steps + 1, functions); i++) {
			table_element_data[i].z = z;
			table_element_data[i].size = functions;
			table_element_data[i].nbase.resize(boost::extents[functions]);
			table_element_data[i].nbase = nbase;
			table_element_data[i].lbase.resize(boost::extents[functions]);
			table_element_data[i].lbase = lbase;
			table_element_data[i].mbase.resize(boost::extents[functions]);
			table_element_data[i].mbase = mbase;
			table_element_data[i].zbase.resize(boost::extents[functions]);
			table_element_data[i].zbase = zbase;
			increment_zbase(functions - 1, steps, min_z, step_z);
			table_element_data[i].coulomb.resize(boost::extents[functions][functions][functions][functions]);
		}
		last_calculation = pow(steps + 1, functions);
		cout << "Starting calculation (" << thread::hardware_concurrency() / concurrency_number << " worker(s))..." << endl;
		for (int i = 0; i < thread::hardware_concurrency() / concurrency_number; i++) {
			worker[i] = thread(thread_runner);
		}
		for (int i = 0; i < thread::hardware_concurrency() / concurrency_number; i++) {
			worker[i].join();
		}
		/* Perform sanity checks on all table points, disregarding clearly faulty results. */
		lowest_energy_sanity_check = -0.5 * pow(z, 3);
		for (int i = 0; i < pow(steps + 1, functions); i++) {
			if (
				(table_element_data[i].total_electron_energy < lowest_energy) &&
				(table_element_data[i].total_electron_energy > lowest_energy_sanity_check) &&
				(abs(table_element_data[i].lowest_eigenvector.maxCoeff()) < 10)
				) {
				lowest_energy = table_element_data[i].total_electron_energy;
				lowest_energy_index = i;
			}
		}
		cout << "Runtime: " << chrono::duration_cast<std::chrono::seconds>(chrono::high_resolution_clock::now() - start_time).count() << " s" << endl;
		cout << "Total electronic energy: " << table_element_data[lowest_energy_index].total_electron_energy << " Ha" << endl;
		cout << "Lowest state effective nuclear charges: ";
		for (int i = 0; i < functions; i++) {
			cout << table_element_data[lowest_energy_index].zbase[i] << " ";
		}
		cout << endl;
		cout << "Lowest state eigenvector: " << endl << table_element_data[lowest_energy_index].lowest_eigenvector << endl;
		cout << "Lowest state eigenvalue: " << table_element_data[lowest_energy_index].lowest_eigenvalue << endl;
		return 0;
	}
	/* User has selected the gradient calculation mode */
	else if (gradient_calculation_mode == true) {
		dataset gradient_data;
		int iterations;
		int z;
		double temp;
		double iteration_size;
		double total_electronic_energy;
		double total_electronic_energy_old = 0;
		boost::multi_array<double, 1> direction;
		cout << "[SETUP] Gradient mode enabled!" << endl;
		cout << "Give the number of basis functions: ";
		cin >> gradient_data.size;
		direction.resize(boost::extents[gradient_data.size]);
		cout << "Give the element number (MUST BE EVEN): ";
		cin >> gradient_data.z;
		gradient_data.nbase.resize(boost::extents[gradient_data.size]);
		gradient_data.lbase.resize(boost::extents[gradient_data.size]);
		gradient_data.mbase.resize(boost::extents[gradient_data.size]);
		gradient_data.zbase.resize(boost::extents[gradient_data.size]); 
		gradient_data.coulomb.resize(boost::extents[gradient_data.size][gradient_data.size][gradient_data.size][gradient_data.size]);
		for (int i = 0; i < gradient_data.size; i++) {
			cout << "Give the N L M quantum numbers of basis function " << i + 1 << ": ";
			cin >> gradient_data.nbase[i];
			cin >> gradient_data.lbase[i];
			cin >> gradient_data.mbase[i];
			cout << "Give the starting point Z of basis function " << i + 1 << ": ";
			cin >> gradient_data.zbase[i];
		}
		cout << "Give the number of iterations: ";
		cin >> iterations;
		cout << "Give the iterative jump length coefficient: ";
		cin >> iteration_size;
		auto start_time = chrono::high_resolution_clock::now();
		/* Calculate the first direction of the negative gradient */
		calculate_gradient(gradient_data, direction);
		for (int i = 0; i < iterations; i++) {
			/* Calculate the total electronic energy for current point */
			total_electronic_energy_old = calculation_runner(gradient_data);
			/* Normalize the direction vector */
			temp = 0;
			for (int j = 0; j < gradient_data.size; j++) {
				temp += pow(direction[j], 2.00);
			}
			temp = pow(temp, -0.5);
			/* Set the direction vector length according to user selection and normalization constant */
			for (int j = 0; j < gradient_data.size; j++) {
				direction[j] = direction[j] * iteration_size * temp;
			}
			/* Make the step in Z-space */
			for (int j = 0; j < gradient_data.size; j++) {
				gradient_data.zbase[j] += direction[j];
			}
			/* Calculate the total electoric energy for current point */
			total_electronic_energy = calculation_runner(gradient_data);
			/* If energy has increased, reevaluate the gradient */
			if (total_electronic_energy > total_electronic_energy_old) {
				calculate_gradient(gradient_data, direction);
			}
		}
		cout << "Runtime: " << chrono::duration_cast<std::chrono::seconds>(chrono::high_resolution_clock::now() - start_time).count() << " s" << endl;
		cout << "Total electronic energy: " << gradient_data.total_electron_energy << " Ha" << endl;
		cout << "Lowest state effective nuclear charges: ";
		for (int i = 0; i < gradient_data.size; i++) {
			cout << gradient_data.zbase[i] << " ";
		}
		cout << endl;
		cout << "Lowest state eigenvector: " << endl << gradient_data.lowest_eigenvector << endl;
		cout << "Lowest state eigenvalue: " << gradient_data.lowest_eigenvalue << endl;
		cout << "Complete linear expansion coefficient matrix: " << endl << gradient_data.all_eigenvectors << endl;
		return 0;
	}
	/* User has selected neither table nor gradient mode. Proceed normally. */
	dataset element_data;
	cout << "Give the amount of basis functions: ";
	cin >> element_data.size;
	cout << "Give the element number (MUST BE EVEN): ";
	cin >> element_data.z;
	element_data.nbase.resize(boost::extents[element_data.size]);
	element_data.lbase.resize(boost::extents[element_data.size]);
	element_data.mbase.resize(boost::extents[element_data.size]);
	element_data.zbase.resize(boost::extents[element_data.size]);
	element_data.coulomb.resize(boost::extents[element_data.size][element_data.size][element_data.size][element_data.size]);
	for (int i = 0; i < element_data.size; i++) {
		cout << "Give the N L M quantum numbers of basis function " << i + 1 << ": ";
		cin >> element_data.nbase[i]; 
		cin >> element_data.lbase[i];
		cin >> element_data.mbase[i];
		cout << "Give the effective nuclear charge of basis function " << i + 1 << ": ";
		cin >> element_data.zbase[i];
	}

	cout << "Total electronic energy: " << calculation_runner(element_data) << " Ha" << endl;
	cout << "Lowest state eigenvector: " << endl << element_data.lowest_eigenvector << endl;
	cout << "Lowest state eigenvalue: " << element_data.lowest_eigenvalue << endl;
	return 0;
}