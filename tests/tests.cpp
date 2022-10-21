#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <murlib/matrix_solver.h>

const double EPSILON = 0.0001;

// Matrix solver algorithms tests
TEST_CASE("the tridiagonal system solver") {
    double A[4] = { 0, 1, 1, 1 };
    double B[4] = { 2, 10, -5, 4 };
    double C[4] = { 1, -5, 2, 0 };
    double D[4] = { -5, -18, -40, -27 };
    double X[4];
    murlib::tridiagonal_matrix(4, A, B, C, D, X);


    double target_x[] = { -3, 1, 5, -8 };
    for (int i = 0; i < 4; ++i) {
        CHECK(fabs(target_x[i] - X[i]) < EPSILON);
    }
}
TEST_CASE("gauss elimination system solver") {
    double A[9] = { 1,	2,	3,
                    2,	2,	1,
                    1,	2,	2 };    
    double b[3] = { 1, 0, 1 };
    double X[3];
    murlib::solve_gauss(3, A, b, X);


    double target_x[] = { -1, 1, 0 };
    for (int i = 0; i < 3; ++i) {
        CHECK(fabs(target_x[i] - X[i]) < EPSILON);
    }
}
TEST_CASE("gauss elimination system solver with zero pivot element") {
    double A[9] = { 0,	2,	1,
                    1,	1,	1,
                    1,	1,	2	};
    double b[3] = { 1, 0, 1 };
    double X[3];
    murlib::solve_gauss(3, A, b, X);


    double target_x[] = { -1, 0, 1 };
    for (int i = 0; i < 3; ++i) {
        CHECK(fabs(target_x[i] - X[i]) < EPSILON);
    }
}
TEST_CASE("gauss elimination system solver when det(A) = 0") {
    double A[9] = { 0,	2,	1,
                    0,	1,	1,
                    0,	1,	2 };
    double b[3] = { 1, 0, 1 };
    double X[3];
    bool is_error = false;
    try {
        murlib::solve_gauss(3, A, b, X);
    }
    catch (std::invalid_argument e) {
        is_error = true;
    }

    CHECK(is_error);

}
TEST_CASE("lu decompostion") {
    double A[9] = { 
        1,	2,	3,
        2,	2,	1,
        1,	2,	2 
    };
    double L[9];
    double U[9];
    murlib::lu_decompostion(3, A, L, U);

    double L_target[9] = {
        1, 0., 0.,
        2, 1, 0.,
        1, 0., 1
    };

    double U_target[9] = {
        1, 2, 3,
        0., -2, -5,
        0., 0., -1
    };

    for (int i = 0; i < 9; ++i) {
        CHECK(abs(L_target[i] - L[i]) < EPSILON);
    }
    for (int i = 0; i < 9; ++i) {
        CHECK(abs(U_target[i] - U[i]) < EPSILON);
    }
}
TEST_CASE("lu system solver") {
    double A[9] = { 1,	2,	3,
                    2,	2,	1,
                    1,	2,	2 };
    double b[3] = { 1, 0, 1 };
    double X[3];

    double L[9];
    double U[9];

    murlib::lu_decompostion(3, A, L, U);
    murlib::solve_lu(3, L, U, b, X);


    double target_x[] = { -1, 1, 0 };
    for (int i = 0; i < 3; ++i) {
        CHECK(fabs(target_x[i] - X[i]) < EPSILON);
    }
}