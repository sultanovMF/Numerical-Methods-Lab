#include <iostream>
#include <vector>
#include <span>
#include<numeric>
#include <effolkronium/random.hpp>

class BandMatrix {
public:
	BandMatrix(size_t N, size_t d) : d(d), N(N) {
		m.resize(N * (d + 1));
		std::fill(m.begin(), m.end(), 0.);
	}

	double operator()(int i, int j) const
	{
		if (i > j) {
			std::swap(i, j);
		}
		if (i < std::max(0, j - d)) {
			return 0.;
		}
		return m[(d + i - j) * N + j];
	}

	int get_d() const {
		return d;
	}

	int get_N() const {
		return N;
	}

	friend void generate_band_matrix(BandMatrix& A, BandMatrix& B, BandMatrix& C);

private:
	std::vector<double> m;
	int d;
	int N;
};

std::ostream& operator<<(std::ostream& os, const BandMatrix& obj)
{
	for (int i = 0; i < obj.get_N(); ++i) {
		for (int j = 0; j < obj.get_N(); ++j) {
			os << std::setprecision(5) << obj(i, j) << " ";
		}
		os << "\n";
	}
	return os;
}

void generate_band_matrix(BandMatrix& A, BandMatrix& B, BandMatrix& C) {
	using Random = effolkronium::random_static;
	for (int i = 0; i < A.d; ++i) {
		for (int j = std::max(0, A.d - i); j < A.N; ++j) {
			double a = Random::get(-1., 1.);
			A.m[i * A.N + j] = a;
			B.m[i * A.N + j] = a;
			C.m[i * A.N + j] = a;
		}
	}

	for (int i = 0; i < A.N; ++i) {
		double sum = 0;
		for (int j = 0; j < A.N; ++j) {
			sum += abs(A(i, j));
		}

		A.m[A.d * A.N + i] = 1.1 * sum;
		B.m[A.d * A.N + i] = 2 * sum;
		C.m[A.d * A.N + i] = 10 * sum;
	}
}

int yacobi_matrix_solver(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, const int m = 10000, double eps = 1e-6) {
	std::vector<double> x_old;
	double error;
	int steps = 0;
	do {
		x_old = x;
		error = 0;
		for (int i = 0; i < N; ++i) {
			x[i] = b[i];
			for (int j = 0; j < i; ++j) {
				x[i] -= A(i, j) * x_old[j];
			}
			for (int j = i + 1; j < N; ++j) {
				x[i] -= A(i, j) * x_old[j];
			}
			x[i] /= A(i, i);

			error = std::max(error, abs(x[i] - x_old[i]));
		}
		++steps;
	} while (error > eps && steps < m);

	return steps;
}

int sor_matrix_solver(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, double w, double eps = 1e-6) {
	std::vector<double> x_old;
	double error;
	int steps = 0;
	do {
		x_old = x;
		error = 0;
		for (int i = 0; i < N; ++i) {
			x[i] = b[i];
			for (int j = 0; j < i; ++j) {
				x[i] -= A(i, j) * x[j];
			}
			for (int j = i + 1; j < N; ++j) {
				x[i] -= A(i, j) * x_old[j];
			}
			x[i] = x[i] * w / A(i, i) + (1 - w) * x_old[i];

			error = std::max(error, abs(x[i] - x_old[i]));
		}
		++steps;
	} while (error > eps);

	return steps;
}
double dot_product(std::vector<double> vect_A, std::vector<double> vect_B)
{
	double product = 0;

	for (int i = 0; i < vect_A.size(); i++)
		product = product + vect_A[i] * vect_B[i];
	return product;
}

int conj_grad_matrix_solver(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, double eps = 1e-6) {
	// init r0
	std::vector<double> r_old(N);
	for (int i = 0; i < N; ++i) {
		r_old[i] = b[i];
		for (int j = 0; j < N; ++j) {
			r_old[i] -= A(i, j) * x[j];
		}
	}
	// init z0
	std::vector<double> z = r_old;
	
	std::vector<double> r(N);
	double alpha;
	double beta;
	double error = 0;
	int steps = 0;
	do {
		error = 0;

		std::vector<double> Az(N);

		for (int j = 0; j < N; ++j) {
			Az[j] = 0;
			for (int k = 0; k < N; ++k) {
				Az[j] += A(j, k) * z[k];
			}
		}
		alpha = dot_product(r_old, r_old) / dot_product(Az, z);

		for (int i = 0; i < N; ++i) {
			double x_old = x[i];
			x[i] = x[i] + alpha * z[i];
			r[i] = r_old[i] - alpha * Az[i];
			error = std::max(error, abs(x[i] - x_old));
		}
		beta = dot_product(r, r) / dot_product(r_old, r_old);


		for (int i = 0; i < N; ++i) {
			z[i] = r[i] + beta * z[i];
		}

		//error = *std::max_element(r.begin(), r.end()) / *std::max_element(b.begin(), b.end());
		//error = eps + 1;
		r_old = r;
		++steps;
	} while (error > eps);
	return steps;
}

int pred_conj_grad_matrix_solver(int N, const BandMatrix& A, const std::vector<double>& b, std::vector<double>& x, int m, double eps = 1e-6) {
	yacobi_matrix_solver(N, A, b, x, m, eps);
	int count = conj_grad_matrix_solver(N, A, b, x, eps);
	return count;
}
int main() {
	using Random = effolkronium::random_static;
	// Индивидуальное задание
	const int N = 750;
	const int l = 18;
	const double eps = 1e-6;
	bool res_print = false;
	int d = l * 2;
	
	// Генерация матриц
	BandMatrix A1(N, d); // q = 1.1
	BandMatrix A2(N, d); // q = 2
	BandMatrix A3(N, d); // q = 10
	generate_band_matrix(A1, A2, A3);

	auto b = Random::get<std::vector>(-5., 5., N);

	// Метод Якоби q = 1.1
	{
		std::vector<double> x(N, 0);
		int steps = yacobi_matrix_solver(N, A1, b, x);
		std::cout << "Yacobi q = 1.1, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Метод Якоби q = 2
	{
		std::vector<double> x(N, 0);
		int steps = yacobi_matrix_solver(N, A2, b, x);
		std::cout << "Yacobi q = 2, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Метод Якоби q = 10
	{
		std::vector<double> x(N, 0);
		int steps = yacobi_matrix_solver(N, A3, b, x);
		std::cout << "Yacobi q = 10, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}

	std::vector<double> omega_arr = {0.1, 0.2, 0.5, 0.7, 0.9, 1.2, 1.4, 1.6, 1.8};
	std::cout << std::endl;
	// Метод SOR q = 1.1
	{
		for (double omega : omega_arr) {
			std::vector<double> x(N, 0);
			int steps = sor_matrix_solver(N, A1, b, x, omega);
			std::cout << "SOR q = 1.1, " << " w = " << omega << " steps = " << steps << std::endl;
		}

	}
	std::cout << std::endl;
	// Метод SOR q = 2
	{
		for (double omega : omega_arr) {
			std::vector<double> x(N, 0);
			int steps = sor_matrix_solver(N, A2, b, x, omega);
			std::cout << "SOR q = 2, " << " w = " << omega << " steps = " << steps << std::endl;
		}
	}
	std::cout << std::endl;
	// Метод SOR q = 10
	{
		for (double omega : omega_arr) {
			std::vector<double> x(N, 0);
			int steps = sor_matrix_solver(N, A3, b, x, omega);
			std::cout << "SOR q = 10, " << " w = " << omega << " steps = " << steps << std::endl;
		}
	}
	std::cout << std::endl;
	// Zeidel q = 1.1
	{
		std::vector<double> x(N, 0);
		int steps = sor_matrix_solver(N, A1, b, x, 1);
		std::cout << "Zeidel q = 1.1, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Zeidel q = 2
	{
		std::vector<double> x(N, 0);
		int steps = sor_matrix_solver(N, A2, b, x, 1);
		std::cout << "Zeidel q = 2, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Zeidel q = 10
	{
		std::vector<double> x(N, 0);
		int steps = sor_matrix_solver(N, A3, b, x, 1);
		std::cout << "Zeidel q = 10, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	// Метод сопряженных градиентов q = 1.1
	{
		std::vector<double> x(N, 0);
		int steps = conj_grad_matrix_solver(N, A1, b, x);
		std::cout << "Conjugate gradient q = 1.1, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Метод сопряженных градиентов q = 2
	{
		std::vector<double> x(N, 0);
		int steps = conj_grad_matrix_solver(N, A2, b, x);
		std::cout << "Conjugate gradient q = 2, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	// Метод сопряженных градиентов q = 10
	{
		std::vector<double> x(N, 0);
		int steps = conj_grad_matrix_solver(N, A3, b, x);
		std::cout << "Conjugate gradient q = 10, steps = " << steps << std::endl;
		if (res_print) {
			std::cout << "Result: ";
			for (int i = 0; i < N; ++i) {
				std::cout << x[i] << " ";
			}

			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;
	// Предобусловленный метод сопряженных градиентов q = 1.1
	{
		for (int m = 5; m < 10; ++m) {
			std::vector<double> x(N, 0);
			int steps = pred_conj_grad_matrix_solver(N, A1, b, x, m);
			std::cout << "Pred conjugate gradient q = 1.1, m = " << m << " steps = " << steps << std::endl;
		}
	}
	// Предобусловленный метод сопряженных градиентов q = 2
	{
		for (int m = 5; m < 10; ++m) {
			std::vector<double> x(N, 0);
			int steps = pred_conj_grad_matrix_solver(N, A2, b, x, m);
			std::cout << "Pred conjugate gradient q = 2, m = " << m << " steps = " << steps << std::endl;
		}
	}
	// Предобусловленный метод сопряженных градиентов q = 10
	{
		for (int m = 5; m < 10; ++m) {
			std::vector<double> x(N, 0);
			int steps = pred_conj_grad_matrix_solver(N, A3, b, x, m);
			std::cout << "Pred conjugate gradient q = 10, m = " << m << " steps = " << steps << std::endl;
		}
	}
}

//std::cout << "\nResult ";
//for (int i = 0; i < N; ++i) {
//	std::cout << x[i] << " ";
//}