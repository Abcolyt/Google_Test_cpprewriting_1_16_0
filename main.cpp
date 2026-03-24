//#define GTEST_HAS_TEST_SUITE_P 1
//#include "pch.h"

// Явные объявления для решения проблемы линковки
namespace testing {
    extern bool FLAGS_gtest_catch_exceptions;
}

bool testing::FLAGS_gtest_catch_exceptions = true;

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include <vector>
#include <utility>
#include <functional>
#include <random>
#include <algorithm>
#include <chrono> 

#if !__has_include("core/complex.h")
#include "file_h/complex.h"
#include "file_h/fraction.h"
#include "file_h/polynomial.h"
#include "file_h/matrix.h"

#include <file_h/interpolators.h>

#include "file_h/integral.h"
#else

#include "core/complex.h"
#include "core/fraction.h"
#include "core/polynomial.h"
#include "linalg/matrix.h"

#include "numerical/.h"

#include "numerical/integration/integral.h"
#endif


/// перенаправление потока вывода
class CoutRedirect {
public:
    CoutRedirect() : old_buf(std::cout.rdbuf()) {
        std::cout.rdbuf(new_buf.rdbuf());
    }

    ~CoutRedirect() {
        std::cout.rdbuf(old_buf);
    }

    // Опционально: получить то, что было "напечатано" во время перенаправления
    std::string getOutput() const {
        return new_buf.str();
    }

private:
    std::streambuf* old_buf;
    std::ostringstream new_buf;
};
/// 
namespace Fraction {
    class FractionTest : public ::testing::Test {
    protected:

    };

    TEST(Constructor, DefaultConstructor) {
        fraction<int> f;
        EXPECT_EQ(f.getNumerator(), 0);
        EXPECT_EQ(f.getDenominator(), 1);
    }

    // Тест конструктора с одним аргументом
    TEST(Constructor, SingleArgConstructor) {
        fraction<int> f(5);
        EXPECT_EQ(f.getNumerator(), 5);
        EXPECT_EQ(f.getDenominator(), 1);
    }

    // Тест конструктора с двумя аргументами
    TEST(Constructor, TwoArgConstructor) {
        fraction<int> f(2, 4);
        EXPECT_EQ(f.getNumerator(), 1);
        EXPECT_EQ(f.getDenominator(), 2);
    }

    // Тест оператора сложения
    TEST(Operator, AdditionOperator) {
        fraction<int> f1(1, 2);
        fraction<int> f2(1, 3);
        fraction<int> result = f1 + f2;
        EXPECT_EQ(result.getNumerator(), 5);
        EXPECT_EQ(result.getDenominator(), 6);
    }

    // Тест оператора умножения
    TEST(Operator, MultiplicationOperator) {
        fraction<int> f1(2, 3);
        fraction<int> f2(3, 4);
        fraction<int> result = f1 * f2;
        EXPECT_EQ(result.getNumerator(), 1);
        EXPECT_EQ(result.getDenominator(), 2);
    }

    // Тест оператора деления
    TEST(Operator, DivisionOperator) {
        fraction<int> f1(2, 3);
        fraction<int> f2(3, 4);
        fraction<int> result = f1 / f2;
        EXPECT_EQ(result.getNumerator(), 8);
        EXPECT_EQ(result.getDenominator(), 9);
    }

    // Тест оператора вычитания
    TEST(Operator, SubtractionOperator) {
        fraction<int> f1(5, 6);
        fraction<int> f2(1, 6);
        fraction<int> result = f1 - f2;
        EXPECT_EQ(result.getNumerator(), 2);
        EXPECT_EQ(result.getDenominator(), 3);
    }

    // Тест метода сокращения дроби
    TEST(Metode, SimplifyFraction) {
        fraction<int> f(8, 12);
        f.reduce();
        EXPECT_EQ(f.getNumerator(), 2);
        EXPECT_EQ(f.getDenominator(), 3);
    }
    // Тест метода сравнения двух дробей
    TEST(Equality_method, CompareFractions) {
        fraction<int> f1(1, 2);
        fraction<int> f2(2, 4);
        bool result = f1 == f2;
        EXPECT_TRUE(result);
    }
}
namespace Polynomial {



    TEST(PolynomialTest, AdditionOperatorTest) {
        polynomial<int> p1; p1.newsize(2);
        p1[0] = 1;
        p1[1] = 2;

        polynomial<int> p2; p2.newsize(2);
        p2[0] = 3;
        p2[1] = 4;

        polynomial<int> sum = p1 + p2;
        EXPECT_EQ(sum[0], 4);
        EXPECT_EQ(sum[1], 6);
    }


    TEST(PolynomialTest, InputOperatorTest) {
        std::stringstream ss("3\n1 2 3\n");
        polynomial<int> p;
        ss >> p;
        EXPECT_EQ(p.get_deg(), 3);
        EXPECT_EQ(p[0], 1);
        EXPECT_EQ(p[1], 2);
        EXPECT_EQ(p[2], 3);
    }

    TEST(OutputOperatorTest, Polynomial_OutputOperatorTest_FULL) {
        polynomial<int> p(3);
        p.newsize(3);
        p.output_mode_set(output_mode::FULL);
        p[0] = 1;
        p[1] = 2;
        p[2] = 3;

        std::stringstream ss;
        ss << p;

        std::string expectedOutput = "Degree: 3, Coefficients: 1 + 2x + 3x^2";
        EXPECT_EQ(ss.str(), expectedOutput);
    }
    TEST(OutputOperatorTest, Polynomial_OutputOperatorTest_ABBREVIATED) {
        polynomial<int> p(3);
        p.newsize(3);
        p.output_mode_set(1);
        p[0] = 1;
        p[1] = 2;
        p[2] = 3;

        std::stringstream ss;
        ss << p;

        std::string expectedOutput = "1+2x+3x^2";
        EXPECT_EQ(ss.str(), expectedOutput);
    }
    TEST(OutputOperatorTest, Polynomial_OutputOperatorTest_SHORT) {
        polynomial<int> p(3);
        p.newsize(3);
        p.output_mode_set(2);
        p[0] = 1;
        p[1] = 2;
        p[2] = 3;

        std::stringstream ss;
        ss << p;

        std::string expectedOutput = "1 2x 3x^2";
        EXPECT_EQ(ss.str(), expectedOutput);
    }

    TEST(PolynomialTest, MultiplyTest) {
        polynomial<int> p1;
        p1.newsize(3);

        p1.output_mode_set(0);
        p1[0] = 2;
        p1[1] = 3;
        p1[2] = 1;

        polynomial<int> p2;
        p2.newsize(3);
        p2[0] = 2;
        p2[1] = 3;
        p2[2] = 1;

        polynomial<int> product; product = p1 * p2;


        std::stringstream ss;
        ss << product;

        std::string expectedOutput1 = "Degree: 5, Coefficients: 4 + 12x + 13x^2 + 6x^3 + 1x^4", expectedOutput2 = "4 12x 13x^2 6x^3 1x^4", expectedOutput3 = "4+12x+13x^2+6x^3+1x^4";
        EXPECT_TRUE(ss.str() == expectedOutput1 || ss.str() == expectedOutput2 || ss.str() == expectedOutput3);


    }

    TEST(PolynomialTest, DivideTest) {
        polynomial<double> p1;
        p1.newsize(3);
        p1[0] = 4;
        p1[1] = 8;
        p1[2] = 12;

        polynomial<double> p2; p2.newsize(1);
        p2[0] = 2;

        polynomial<double> quotient = p1 / p2;
        EXPECT_EQ(quotient.get_deg(), 3);
        EXPECT_EQ(quotient[0], 2);
        EXPECT_EQ(quotient[1], 4);
        EXPECT_EQ(quotient[2], 6);
    }

    TEST(PolynomialTest, ModulusTest) {
        polynomial<double> p1;
        p1.newsize(4);
        p1[0] = 5;
        p1[1] = 10;
        p1[2] = 15;
        p1[3] = 20;

        polynomial<double> p2;
        p2.newsize(2);
        p2[0] = 3;
        p2[1] = 6;

        polynomial<double> remainder = p1 % p2;
        EXPECT_EQ(remainder.get_deg(), 1);
        EXPECT_EQ(remainder[0], 1.25);

    }

}
namespace Matrix {
    TEST(Constructor, ConstructorWithSize)
    {
        matrix<int> m(3, 3);
        EXPECT_EQ(m.getcol(), 3);
        EXPECT_EQ(m.getrow(), 3);
    }


    TEST(operators, MatrixMultiplication_operator1)
    {
        matrix<int> m1(3, 3);
        m1[0][0] = 1; m1[0][1] = 2; m1[0][2] = 3;
        m1[1][0] = 4; m1[1][1] = 5; m1[1][2] = 6;
        m1[2][0] = 7; m1[2][1] = 8; m1[2][2] = 9;

        matrix<int> m2(3, 3);
        m2[0][0] = 9; m2[0][1] = 8; m2[0][2] = 7;
        m2[1][0] = 6; m2[1][1] = 5; m2[1][2] = 4;
        m2[2][0] = 3; m2[2][1] = 2; m2[2][2] = 1;

        matrix<int> result = m1 * m2;

        EXPECT_EQ(result[0][0], 30); EXPECT_EQ(result[0][1], 24); EXPECT_EQ(result[0][2], 18);
        EXPECT_EQ(result[1][0], 84); EXPECT_EQ(result[1][1], 69); EXPECT_EQ(result[1][2], 54);
        EXPECT_EQ(result[2][0], 138); EXPECT_EQ(result[2][1], 114); EXPECT_EQ(result[2][2], 90);
    }

    TEST(operators, MatrixMultiplication_operator2)
    {
        matrix<int> m1(3, 3);
        m1[0][0] = 1; m1[0][1] = 2; m1[0][2] = 3;
        m1[1][0] = 4; m1[1][1] = 5; m1[1][2] = 6;
        m1[2][0] = 7; m1[2][1] = 8; m1[2][2] = 9;

        matrix<int> m2(3, 2);
        m2[0][0] = 9; m2[0][1] = 8;
        m2[1][0] = 6; m2[1][1] = 5;
        m2[2][0] = 3; m2[2][1] = 2;

        matrix<int> result = m1 * m2;
        std::stringstream method_ans;
        method_ans << result;

        //EXPECT_EQ(method_ans.str(), "sizex:3sizey:2\n[0][0] = \t30\t | [0][1] = \t24\t | \n[1][0] = \t84\t | [1][1] = \t69\t | \n[2][0] = \t138\t | [2][1] = \t114\t | \n");

        EXPECT_EQ(result[0][0], 30); EXPECT_EQ(result[0][1], 24);
        EXPECT_EQ(result[1][0], 84); EXPECT_EQ(result[1][1], 69);
        EXPECT_EQ(result[2][0], 138); EXPECT_EQ(result[2][1], 114);
    }



    TEST(methods, MatrixTranspose)//there were no problems in the rest of the program
    {
        matrix<double> m(2, 3);
        m[0][0] = 1; m[0][1] = 2; m[0][2] = 3;
        m[1][0] = 4; m[1][1] = 5; m[1][2] = 6;
        //std::cout << m;
        matrix<double> transposed = m.transpose();

        //std::cout << transposed;
        EXPECT_EQ(transposed.getcol(), 3);
        EXPECT_EQ(transposed.getrow(), 2);
        EXPECT_EQ(transposed[0][0], 1);
        EXPECT_EQ(transposed[1][0], 2);
        EXPECT_EQ(transposed[2][0], 3);
        EXPECT_EQ(transposed[0][1], 4);
        EXPECT_EQ(transposed[1][1], 5);
        EXPECT_EQ(transposed[2][1], 6);
    }

    TEST(methods, determinant) {
        CoutRedirect redirect;
        std::stringstream ss("2\n2\n2\n1\n1\n1\n1\n2\n1\n2\n1\n1\n2\n1\n3\n1\n1\n4\n1\n2\n3\n4\n1\n1\n");
        matrix<fraction<polynomial<int>>> mtrx;
        ss >> mtrx;
        std::stringstream method_ans;
        method_ans << mtrx.determinant();
        std::string true_ans_str1 = "(Degree: 5, Coefficients: 0 + (-2)x + (-1)x^2 + 7x^3 + 4x^4) / (Degree: 1, Coefficients: 1)",
            true_ans_str2 = "(0+(-2)x+(-1)x^2+7x^3+4x^4) / (1)",
            true_ans_str3 = "(0 -2x -1x^2 7x^3 4x^4) / (1)";
        std::cout << mtrx.determinant();
        EXPECT_TRUE(method_ans.str() == true_ans_str1 || method_ans.str() == true_ans_str2 || method_ans.str() == true_ans_str3);
    }

    TEST(methods, inverse_M_int) {
        CoutRedirect redirect;
        std::stringstream ss("2\n2\n2\n1\n1\n1\n1\n2\n1\n2\n1\n1\n2\n1\n3\n1\n1\n4\n1\n2\n3\n4\n1\n1\n");
        matrix<fraction<polynomial<int>>> mtrx;//there used to be a problem with %, and it only worked with int, now it seems to have disappeared
        ss >> mtrx;
        std::stringstream method_ans;
        method_ans << mtrx.inverse_M();
        //"cols: 2, rows: 2\n                                                                  (1 4x 8x^2 12x^3 11x^4 4x^5) / (0 -2x -5x^2 3x^3 17x^4 15x^5 4x^6)                                              (-1 -3x -2x^2) / (0 -2x -3x^2 6x^3 11x^4 4x^5)\n                                                                  (-1 -4x -3x^2) / (0 -2x -3x^2 6x^3 11x^4 4x^5)                                              (1 1x) / (0 -2x -1x^2 7x^3 4x^4)\n"

        std::string true_ans_str1 = "sizex:2sizey:2\n[0][0] = \t(Degree: 6, Coefficients: 1 + 4x + 8x^2 + 12x^3 + 11x^4 + 4x^5) / (Degree: 7, Coefficients: 0 + (-2)x + (-5)x^2 + 3x^3 + 17x^4 + 15x^5 + 4x^6)\t | [0][1] = \tDegree: 3, Coefficients: (-1) + (-3)x + (-2)x^2 / (Degree: 6, Coefficients: 0 + (-2)x + (-3)x^2 + 6x^3 + 11x^4 + 4x^5)\t | \n[1][0] = \tDegree: 3, Coefficients: (-1) + (-4)x + (-3)x^2 / (Degree: 6, Coefficients: 0 + (-2)x + (-3)x^2 + 6x^3 + 11x^4 + 4x^5)\t | [1][1] = \tDegree: 2, Coefficients: 1 + 1x / (Degree: 5, Coefficients: 0 + (-2)x + (-1)x^2 + 7x^3 + 4x^4)\t | \n",
            true_ans_str2 = "cols: 2, rows: 2\n                                                                  (1 4x 8x^2 12x^3 11x^4 4x^5) / (0 -2x -5x^2 3x^3 17x^4 15x^5 4x^6)                                              (-1 -3x -2x^2) / (0 -2x -3x^2 6x^3 11x^4 4x^5)\n                                                                  (-1 -4x -3x^2) / (0 -2x -3x^2 6x^3 11x^4 4x^5)                                              (1 1x) / (0 -2x -1x^2 7x^3 4x^4)\n",
            true_ans_str3 = "sizex:2sizey:2\n[0][0] = \t1 4x 8x^2 12x^3 11x^4 4x^5 / (0 -2x -5x^2 3x^3 17x^4 15x^5 4x^6)\t | [0][1] = \t-1 -3x -2x^2 / (0 -2x -3x^2 6x^3 11x^4 4x^5)\t | \n[1][0] = \t-1 -4x -3x^2 / (0 -2x -3x^2 6x^3 11x^4 4x^5)\t | [1][1] = \t1 1x / (0 -2x -1x^2 7x^3 4x^4)\t | \n";
        EXPECT_EQ(method_ans.str(), true_ans_str2);
        EXPECT_TRUE(method_ans.str() == true_ans_str1 || method_ans.str() == true_ans_str2 || method_ans.str() == true_ans_str3);
    }
    TEST(methods, inverse_M_float) {
        CoutRedirect redirect;
        std::stringstream ss("2\n2\n2\n1\n1\n1\n1\n2\n1\n2\n1\n1\n2\n1\n3\n1\n1\n4\n1\n2\n3\n4\n1\n1\n");
        matrix<fraction<polynomial<float>>> mtrx;//there used to be a problem with %, and it only worked with int, now it seems to have disappeared
        ss >> mtrx;
        std::stringstream method_ans;
        method_ans << mtrx.inverse_M();
        std::string true_ans_str1 = "cols: 2, rows: 2\n                                                                                                                                                               (-233827 -935306x -1.87061e+06x^2 -2.80592e+06x^3 -2.57209e+06x^4 -935306x^5) / (0 467653x 1.16913e+06x^2 -701480x^3 -3.97505e+06x^4 -3.5074e+06x^5 -935306x^6)                                                                                       (36.125 108.375x 72.2499x^2) / (0 72.2499x 108.375x^2 -216.75x^3 -397.375x^4 -144.5x^5)\n                                                                                                                                                               (-2.89286 -8.67857x) / (-0 -5.78571x -2.89286x^2 20.25x^3 11.5714x^4)                                                                                       (-0.5 -0.5x) / (0 1x 0.5x^2 -3.5x^3 -2x^4)\n";

        EXPECT_EQ(method_ans.str(), true_ans_str1);
    }
    TEST(methods, inverse_M_double) {
        CoutRedirect redirect;
        std::stringstream ss("2\n2\n2\n1\n1\n1\n1\n2\n1\n2\n1\n1\n2\n1\n3\n1\n1\n4\n1\n2\n3\n4\n1\n1\n");
        matrix<fraction<polynomial<double>>> mtrx;//there used to be a problem with %, and it only worked with int, now it seems to have disappeared
        ss >> mtrx;
        std::stringstream method_ans;
        method_ans << mtrx.inverse_M();

        std::string true_ans_str1 = "cols: 2, rows: 2\n                                                                                                                               (-36.125 -144.5x -289x^2 -433.5x^3 -397.375x^4 -144.5x^5) / (0 72.25x 180.625x^2 -108.375x^3 -614.125x^4 -541.875x^5 -144.5x^6)                                                                                                                             (-2.82452e+14 -8.47357e+14x -5.64905e+14x^2) / (0 -5.64905e+14x -8.47357e+14x^2 1.69471e+15x^3 3.10698e+15x^4 1.12981e+15x^5)\n                                                                                                                               (-4.5036e+15 -1.80144e+16x -1.35108e+16x^2) / (-0 -9.0072e+15x -1.35108e+16x^2 2.70216e+16x^3 4.95396e+16x^4 1.80144e+16x^5)                                                                                                                             (-0.5 -0.5x) / (0 1x 0.5x^2 -3.5x^3 -2x^4)\n";
        std::cout << method_ans.str();
        EXPECT_EQ(method_ans.str(), true_ans_str1);
    }


    TEST(EigenvaluesTest, simple_QR_decomposition) {
        const int n = 4;
        const double eps = 1e-5;
        matrix<double> A_diag = matrix<double>::randomDiagonal(n, -100.0, 100.0);

        std::vector<double> diag_values;
        for (int i = 0; i < n; ++i) {
            diag_values.push_back(A_diag[i][i]);
        }
        std::sort(diag_values.begin(), diag_values.end());
        matrix<double> T, T_inv;
        bool is_invertible = false;
        int attempts = 0;

        while (!is_invertible && attempts < 10) {
            T = matrix<double>::random(n, n, -100.0, 100.0);
            try {
                T_inv = T.inverse_M();
                is_invertible = true;
            }
            catch (...) {
                attempts++;
            }
        }

        if (!is_invertible) {
            FAIL() << "Не удалось сгенерировать обратиую матрицу T";
        }

        matrix<double> A = T * A_diag * T_inv;

        // Вычисление собственных значений
        std::vector<std::complex<double>> eigenvalues = matrixfunction::ComputeEigenvalues(A);
        ASSERT_EQ(eigenvalues.size(), n) << "Количество собственных значений не совпадает";

        // Проверка мнимых частей и сбор вещественных частей
        std::vector<double> real_parts;
        for (const auto& ev : eigenvalues) {
            ASSERT_NEAR(ev.imag(), 0.0, eps) << "Обнаружена значимая мнимая часть";
            real_parts.push_back(ev.real());
        }

        std::sort(real_parts.begin(), real_parts.end());

        // Сравнение с исходными значениями
        for (int i = 0; i < n; ++i) {
            ASSERT_NEAR(real_parts[i], diag_values[i], eps)
                << "Отличие на позиции " << i << ": " << real_parts[i] << " vs " << diag_values[i];
        }

    }

#if 0
    namespace ParamEigenvaluesTests {



        using EigenFunc = std::function<std::vector<std::complex<double>>(const matrix<double>&)>;

        struct TestParams {
            int matrix_size;
            EigenFunc eigen_func;
            std::string func_name;
        };

        class EigenvaluesTest : public testing::TestWithParam<TestParams> {
        protected:
            void SetUp() override {
                params = GetParam();
            }

            void runTest() {
                const int n = params.matrix_size;
                const double eps = 1e-5;

                matrix<double> A_diag = matrix<double>::randomDiagonal(n, -100.0, 100.0);

                std::vector<double> diag_values;
                for (int i = 0; i < n; ++i) {
                    diag_values.push_back(A_diag[i][i]);
                }
                std::sort(diag_values.begin(), diag_values.end());

                matrix<double> T, T_inv;
                bool is_invertible = false;
                int attempts = 0;
                const int max_attempts = 10;

                while (!is_invertible && attempts < max_attempts) {
                    T = matrix<double>::random(n, n, -100.0, 100.0);
                    try {
                        T_inv = T.inverse_M();
                        is_invertible = true;
                    }
                    catch (...) {
                        attempts++;
                    }
                }

                if (!is_invertible) {
                    FAIL() << "Failed to generate invertible matrix T after "
                        << max_attempts << " attempts";
                }

                matrix<double> A = T * A_diag * T_inv;

                std::vector<std::complex<double>> eigenvalues = params.eigen_func(A);
                ASSERT_EQ(eigenvalues.size(), n)
                    << "Number of eigenvalues doesn't match matrix size";



                std::vector<double> real_parts, image_part;

                for (const auto& ev : eigenvalues) {
                    ASSERT_NEAR(ev.imag(), 0.0, eps)
                        << "Significant imaginary part found in eigenvalue";
                    real_parts.push_back(ev.real());
                    //image_part.
                }


                std::sort(real_parts.begin(), real_parts.end());
                for (int i = 0; i < n; ++i) {
                    EXPECT_NEAR(real_parts[i], diag_values[i], eps)
                        << "Eigenvalue mismatch at position " << i
                        << " for function " << params.func_name
                        << " and matrix size " << n;
                }


                std::sort(real_parts.begin(), real_parts.end());
                for (int i = 0; i < n; ++i) {
                    EXPECT_NEAR(real_parts[i], diag_values[i], eps)
                        << "Eigenvalue mismatch at position " << i
                        << " for function " << params.func_name
                        << " and matrix size " << n;
                }
            }

            TestParams params;
        };




        std::vector<std::complex<double>> wrap_compute_eigenvalues(const matrix<double>& m) {
            return matrixfunction::ComputeEigenvalues(m);
        }


        std::ostream& operator<<(std::ostream& os, const TestParams& params) {
            os << params.func_name << "_Size_" << params.matrix_size;
            return os;
        }

        struct TestParamNameGenerator {
            template <class ParamType>
            std::string operator()(const testing::TestParamInfo<ParamType>& info) {
                std::ostringstream oss;
                oss << info.param;
                return oss.str();
            }
        };

        TEST_P(EigenvaluesTest, SimilarityTransformation) {
            runTest();
        }


        INSTANTIATE_TEST_CASE_P(
            Matrix_EigenTests,
            EigenvaluesTest,
            ::testing::Values(
                TestParams{ 4 , wrap_compute_eigenvalues, "eig" },
                TestParams{ 7 , [](auto& m) { return matrixfunction::ComputeEigenvalues(m); }, "eig" },
                TestParams{ 50, [](auto& m) { return matrixfunction::ComputeEigenvalues(m); }, "eig" },
                TestParams{ 10, [](auto& m) { return matrixfunction::ComputeEigenvalues3QR(m); }, "eig_3_qr" },
                TestParams{ 50, [](auto& m) { return matrixfunction::ComputeEigenvalues3QR(m); }, "eig_3_qr" }
            ),
            TestParamNameGenerator()
        );
    }
#else
    namespace ParamEigenvaluesTests {

        using EigenFunc = std::function<std::vector<std::complex<double>>(const matrix<double>&)>;

        struct TestParams {
            int matrix_size;
            EigenFunc eigen_func;
            std::string func_name;
        };

        class EigenvaluesTest : public testing::TestWithParam<TestParams> {
        protected:
            void SetUp() override {
                params = GetParam();
            }

            void runTest() {
                const int n = params.matrix_size;
                const double eps = 1e-5;

                matrix<double> A_diag = matrix<double>::randomDiagonal(n, -100.0, 100.0);

                std::vector<double> diag_values;
                for (int i = 0; i < n; ++i) {
                    diag_values.push_back(A_diag[i][i]);
                }
                std::sort(diag_values.begin(), diag_values.end());

                matrix<double> T, T_inv;
                bool is_invertible = false;
                int attempts = 0;
                const int max_attempts = 10;

                while (!is_invertible && attempts < max_attempts) {
                    T = matrix<double>::random(n, n, -100.0, 100.0);
                    try {
                        T_inv = T.inverse_M();
                        is_invertible = true;
                    }
                    catch (...) {
                        attempts++;
                    }
                }

                if (!is_invertible) {
                    FAIL() << "Failed to generate invertible matrix T after "
                        << max_attempts << " attempts";
                }

                matrix<double> A = T * A_diag * T_inv;

                std::vector<std::complex<double>> eigenvalues = params.eigen_func(A);
                ASSERT_EQ(eigenvalues.size(), n)
                    << "Number of eigenvalues doesn't match matrix size";

                std::vector<double> real_parts, image_part;

                for (const auto& ev : eigenvalues) {
                    ASSERT_NEAR(ev.imag(), 0.0, eps)
                        << "Significant imaginary part found in eigenvalue";
                    real_parts.push_back(ev.real());
                }

                std::sort(real_parts.begin(), real_parts.end());
                for (int i = 0; i < n; ++i) {
                    EXPECT_NEAR(real_parts[i], diag_values[i], eps)
                        << "Eigenvalue mismatch at position " << i
                        << " for function " << params.func_name
                        << " and matrix size " << n;
                }
            }

            TestParams params;
        };

        std::vector<std::complex<double>> wrap_compute_eigenvalues(const matrix<double>& m) {
            return matrixfunction::ComputeEigenvalues(m);
        }

        std::ostream& operator<<(std::ostream& os, const TestParams& params) {
            os << params.func_name << "_Size_" << params.matrix_size;
            return os;
        }

        struct TestParamNameGenerator {
            template <class ParamType>
            std::string operator()(const testing::TestParamInfo<ParamType>& info) {
                std::ostringstream oss;
                oss << info.param;
                return oss.str();
            }
        };

        TEST_P(EigenvaluesTest, SimilarityTransformation) {
            runTest();
        }

        // заменено на INSTANTIATE_TEST_SUITE_P
        INSTANTIATE_TEST_SUITE_P(
            Matrix_EigenTests,
            EigenvaluesTest,
            ::testing::Values(
                TestParams{ 4 , wrap_compute_eigenvalues, "eig" },
                TestParams{ 7 , [](auto& m) { return matrixfunction::ComputeEigenvalues(m); }, "eig" },
                TestParams{ 50, [](auto& m) { return matrixfunction::ComputeEigenvalues(m); }, "eig" },
                TestParams{ 10, [](auto& m) { return matrixfunction::ComputeEigenvalues3QR(m); }, "eig_3_qr" },
                TestParams{ 50, [](auto& m) { return matrixfunction::ComputeEigenvalues3QR(m); }, "eig_3_qr" }
            ),
            TestParamNameGenerator()
        );
    }
#endif


}

//
//#include <boost/safe_numerics/safe_integer.hpp>

namespace Polynomial_counting_methods_2 {
    TEST(nuton, first_static_data1) {
        using namespace polynomial_interpolation;
        std::vector<std::pair<int, int>> Array_xy = { {1, 10}, {2, 20}, {3, 30},  {4,40},{5,50},{6,60},{7,70},{8,80},
            {1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30},{1, 10}, {2, 20}, {3, 30} };

        std::stringstream local_ans, true_ans;
        local_ans << NutonInterpolation(Array_xy);

        true_ans << "0 10x";

        EXPECT_EQ(local_ans.str(), true_ans.str());
    }
    TEST(nuton, second_static_data) {
        using namespace polynomial_interpolation;
        std::vector<std::pair<int, int>> Array_xy = {
        {1, 11}, {2, 21}, {3, 31}, {4,41}, {5,51},
        {1, 999}, {2, 888}, {3, 777},
        {1, 1000}, {2, 2000}, {3, 3000}
        };
        std::stringstream local_ans, true_ans;
        local_ans << NutonInterpolation(Array_xy);

        true_ans << "1 10x";

        EXPECT_EQ(local_ans.str(), true_ans.str());
    }

    // Generator of the vector x y via a lambda function
    template<typename T, typename Func>    std::vector<std::pair<T, T>> generatePointsLambda(int k, T x0, T step, Func F) {
        std::vector<std::pair<T, T>> points;
        for (int i = 0; i < k; ++i) {
            T x = x0 + i * step;
            points.emplace_back(x, F(x));
        }
        return points;
    }
    TEST(nuton, third_static_data) {
        using namespace polynomial_interpolation;


        auto Func = [](float x) { return  1 + 10 * x + 10 * x * x + 10 * x * x * x + 10 * x * x * x * x; };
        auto Array_xy = generatePointsLambda(6, -4, 1, Func);

        std::stringstream local_ans, true_ans;
        local_ans << NutonInterpolation(Array_xy);

        true_ans << "1 10x 10x^2 10x^3 10x^4";

        EXPECT_EQ(local_ans.str(), true_ans.str());
    }

    using namespace polynomial_interpolation;
    using namespace technical_functions;
    using namespace point_generators;
    TEST(nuton, first_dinamic_data_int) {
        polynomial<int> pol;
        pol = GenerateRandomIntCoefficients(3, 9, -20, 20);
        //std::cout << pol << '\n';
        auto Array_xy = GeneratePointsFuncPtr<int>(pol.get_deg() + 3, -4, 1, pol, std::function<int(polynomial<int>, int)>(polynomialfunctions::f_polyn_x0_<int>));

        using namespace polynomial_interpolation;
        //std::cout << NutonInterpolation(Array_xy);


        std::stringstream local_ans, true_ans;
        local_ans << NutonInterpolation(Array_xy);
        true_ans << pol;
        EXPECT_EQ(local_ans.str(), true_ans.str());
    }
    TEST(nuton, Array_dinamic_data_int) {
        for (size_t i = 0; i < 25; i++)
        {
            using Type = int;
            polynomial<Type> pol;
            pol = GenerateRandomIntCoefficients(3, 8, -10, 10);
            //std::cout << pol << '\n';
            auto Array_xy = GeneratePointsFuncPtr<Type>(pol.get_deg() + 3, -4, 1, pol, std::function<Type(polynomial<Type>, Type)>(polynomialfunctions::f_polyn_x0_<Type>));

            using namespace polynomial_interpolation;
            //std::cout << NutonInterpolation(Array_xy);


            std::stringstream local_ans, true_ans;
            local_ans << NutonInterpolation(Array_xy);
            pol = pol.cutbag();
            true_ans << pol;
            EXPECT_EQ(local_ans.str(), true_ans.str());
        }

    }

}



//1) тип данных интегриования
//2) тип формулы ньютона котеса

//2.1) интегрируемая функция
//2.2) интервал интегрирования
//2.3) число разбиений интервала на каждом из который будет происходить применение соответсвующей формулы.


// Перечисление для методов Ньютона-Котеса
namespace Integral_Computing_Function {

    //using namespace counting_methods_3;

    // Базовая структура для тестовых данных 
    template<typename TArg, typename TResult>
    struct TestConfig {
        std::function<TResult(TArg)> function;      // Интегрируемая функция
        std::function<TResult(TArg)> reference;     // Первообразная (эталонная функция)
        TArg a;
        TArg b;
        TResult expected;                           // Ожидаемый результат
        TResult tolerance;                          // Допустимая погрешность
        int points;
        std::string name;                           // Имя конфигурации( для отладки)
        TArg alpha;                                 // Дополнительный параметр
        TArg beta;
        int max_iteration_multiplication = 10;
    };

    //конфиги
    namespace configuration {
        // Тестовые конфигурации для double
        std::vector<TestConfig<double, double>> double_configs = {
            {
                [](double x) { return x * x; },           // f(x) = x²
                [](double x) { return x * x * x / 3.0; }, // F(x) = x³/3
                0.0, 1.0,                                 // [0, 1]
                1.0 / 3.0,                                  // ожидаемый ∫x²dx = 1/3
                1e-9,                                     // допуск
                100,                                      // начальное количество точек
                "x_squared_0_1",
                0,0,40
            },
            {
                [](double x) { return std::sin(x); },     // f(x) = sin(x)
                [](double x) { return -std::cos(x); },    // F(x) = -cos(x)
                0.0, M_PI,                                // [0, π]
                2.0,                                      // ожидаемый ∫sin(x)dx = 2
                1e-9,
                100,
                "sin_0_pi",
                0,0,20
            },
            {
                [](double x) { return 1.0 / x; },         // f(x) = 1/x
                [](double x) { return std::log(x); },     // F(x) = ln(x)
                1.0, M_E,                                 // [1, e]
                1.0,                                      // ожидаемый ∫(1/x)dx = 1
                1e-9,
                200,
                "inverse_x_1_to_e",
                0,0,40
            }
        };

        // Тестовые конфигурации для float  
        std::vector<TestConfig<float, float>> float_configs = {
            {
                [](float x) { return x; },
                [](float x) { return x * x / 2.0f; },
                0.0f, 1.0f,
                0.5f,
                1e-3f,
                50,
                "linear_float",
                0,0,10
            }
        };
        using namespace std;
         
       
       //integral =  ((exp(2)) ^ (1 / 5)* (103296 * sin(144) - 37186560 * cos(144)) + (e ^ 2) ^ (1 / 15) * (103554360 * e ^ 2 * sin(112) + 1972464 * e ^ 2 * cos(112)) + (e ^ 7) ^ (1 / 80) * (37186560 * cos(63 / 2) - 103296 * sin(63 / 2)) + (e ^ 7) ^ (1 / 15) * (-103554360 * sin(49 / 2) - 1972464 * cos(49 / 2)) + 6798220455) / (278901352)
        
       std::vector<TestConfig<double, double>> variant_configs = {
            {
                [](double x) { return x * x; },           // f(x) = x²
                [](double x) { return x * x * x / 3.0; }, // F(x) = x³/3
                0.0, 1.0,                                 // [0, 1]
                1.0 / 3.0,                                  // ожидаемый ∫x²dx = 1/3
                1e-9,                                     // допуск
                100,                                      // начальное количество точек
                "x_squared_0_1",
                0,0,
                1
            },
            {
                [](double x) { return std::sin(x); },     // f(x) = sin(x)
                [](double x) { return -std::cos(x); },    // F(x) = -cos(x)
                0.0, M_PI,                                // [0, π]
                2.0,                                      // ожидаемый ∫sin(x)dx = 2
                1e-9,
                100,
                "sin_0_pi",
                0,0,
                1
            },
            {    
               [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
               [](double x) { return (-109680*sin((9*x)/2)-3948480*cos((9*x)/2)+1062243*exp((19*x)/(24))*sin((7*x)/2)+202332*exp((19*x)/(24))*cos((7*x)/2))/(2963645*exp(x/8))+(5*x*x)/2+3746148/(2963645) ; },
               0.7, 3.2,
               20.235345,
               1e-9,
               30,
               "variant_10_alpha_beta_0",
               0.0,0.0,
               25
            },
            {
               [](double x) { return 1.3 * cos(3.5 * x) * exp(2 * x / 3) + 6 * sin(4.5 * x) * exp(-x / 8) + 5 * x;  },
               [](double x) { return (-109680 * sin((9 * x) / 2) - 3948480 * cos((9 * x) / 2) + 1062243 * exp((19 * x) / (24)) * sin((7 * x) / 2) + 202332 * exp((19 * x) / (24)) * cos((7 * x) / 2)) / (2963645 * exp(x / 8)) + (5 * x * x) / 2 + 3746148 / (2963645); },
               0.7, 3.2,
               24.142092678433,
               1e-9,
               30,
               "variant_10_alpha_0_beta_0_and_75",
               0.0, 1.0/4,
               3
            },
            {
               [](double x) { return 1.5* cos(3.7*x) *exp(4*x/ 7) + 3 *sin(2.5*x)* exp(3*x/ 4) + 3*x;  },
               [](double x) { return (27195*exp((4*x)/7)*sin((3.7*x)))/(68681)+(4200*exp((4*x)/7)*cos(3.7*x))/(68681)+(36*exp(0.75*x)*sin(2.5*x))/(109)-(120*exp(0.75*x)*cos(2.5*x))/(109)+(3*x*x)/2+7783920/(7486229); },
               1.5, 3.0,
               5.6085702,
               1e-9,
               30,
               "variant_20_alpha_beta_0",
               0.0,0.0,
               3
            },
            {
               [](double x) { return 1.5 * cos(3.7 * x) * exp(4 * x / 7) + 3 * sin(2.5 * x) * exp(3 * x / 4) + 3 * x;  },
               [](double x) { return (27195 * exp((4 * x) / 7) * sin((3.7 * x))) / (68681) + (4200 * exp((4 * x) / 7) * cos(3.7 * x)) / (68681) + (36 * exp(0.75 * x) * sin(2.5 * x)) / (109) - (120 * exp(0.75 * x) * cos(2.5 * x)) / (109) + (3 * x * x) / 2 + 7783920 / (7486229); },
               1.5, 3.0,
               5.6085702,
               1e-9,
               30,
               "variant_20_alpha_0_beta_5__6",
               0,5.0/6,
               3
            }
        };
       
    }
    using namespace configuration;

    template<typename TArg, typename TResult>
    TResult compute_expected(const TestConfig<TArg, TResult>& config) {
        return config.reference(config.b) - config.reference(config.a);
    }

    template<typename TResult>
    bool check_accuracy(TResult result, TResult expected, TResult tolerance) {
        if constexpr (std::is_integral_v<TResult>) {
            return result == expected;
        }
        else if constexpr (std::is_same_v<TResult, std::complex<double>>) {
            return (std::abs(result.real() - expected.real()) <= tolerance.real()) &&
                (std::abs(result.imag() - expected.imag()) <= tolerance.imag());
        }
        else {
            return std::abs(result - expected) <= tolerance;
        }
    }

    // Шаблонный класс для параметризованных тестов с методом как шаблонным параметром
    template<typename TArg, typename TResult, IntegrateMethod Method>
    class NewtonCotesParamTest : public ::testing::TestWithParam<TestConfig<TArg, TResult>> {
    protected:
        void RunTest() {
            const auto& config = this->GetParam();
            TResult expected = compute_expected(config);

            int current_points = config.points;
            TResult result;
            bool accuracy_achieved = false;
            int iterations = 0;
            //const int max_iterations = 10;

            for (; iterations < config.max_iteration_multiplication; ++iterations) {
                result = unsafe_integrate<Method, TResult, TArg>(
                    config.function, config.a, config.b, current_points
                );
                
                if (check_accuracy(result, expected, config.tolerance)) {
                    accuracy_achieved = true;
                    break;
                }

                current_points *= 2;
            }

            if constexpr (std::is_integral_v<TResult>) {
                EXPECT_EQ(result, expected);
            }
            else if constexpr (std::is_same_v<TResult, std::complex<double>>) {
                EXPECT_NEAR(result.real(), expected.real(), config.tolerance.real());
                EXPECT_NEAR(result.imag(), expected.imag(), config.tolerance.imag());
            }
            else {
                EXPECT_NEAR(result, expected, config.tolerance);
            }
        }
    };

    // Макрос для регистрации тестов для конкретного метода
#define REGISTER_INTEGRATE_TEST(TestName, Method, ArgType, ResultType, ...) \
        using TestName##_##Method##_##ArgType##_##ResultType = \
            NewtonCotesParamTest<ArgType, ResultType, IntegrateMethod::Method>; \
        TEST_P(TestName##_##Method##_##ArgType##_##ResultType, TestName) { \
            RunTest(); \
        } \
        INSTANTIATE_TEST_SUITE_P(TestName, \
            TestName##_##Method##_##ArgType##_##ResultType, \
            ::testing::Values(__VA_ARGS__), \
            [](const ::testing::TestParamInfo<TestName##_##Method##_##ArgType##_##ResultType::ParamType>& info) { \
                return info.param.name + "_" + #Method + "_initial_points_" + std::to_string(info.param.points); \
            })

#define REGISTER_ALL_COTES_INTEGRATE_METHODS_FOR_TYPE(TestName, ArgType, ResultType, ...) \
        REGISTER_INTEGRATE_TEST(TestName, LEFT_RECTANGLE, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, MIDDLE_RECTANGLE, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, TRAPEZOID, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, SIMPSON, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_4_POINT , ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_5_POINT , ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_6_POINT , ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_7_POINT , ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_8_POINT , ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, NEWTON_COTES_9_POINT , ArgType, ResultType, __VA_ARGS__) 

#define REGISTER_ALL_GAUS_INTEGRATE_METHODS_FOR_TYPE(TestName, ArgType, ResultType, ...) \
        REGISTER_INTEGRATE_TEST(TestName, GAUS_3_POINT, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, GAUS_4_POINT, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, GAUS_5_POINT, ArgType, ResultType, __VA_ARGS__); \
        REGISTER_INTEGRATE_TEST(TestName, GAUS_6_POINT, ArgType, ResultType, __VA_ARGS__)

    REGISTER_INTEGRATE_TEST(SpecificSimpsonTest, SIMPSON, double, double,
        configuration::double_configs[0]  // Только x² для метода Симпсона
    );

    // Регистрация всех методов для double
    REGISTER_ALL_COTES_INTEGRATE_METHODS_FOR_TYPE(DoubleTests, double, double,
        double_configs[0], double_configs[1], double_configs[2], variant_configs[2]
    );

    // Регистрация всех методов для float
    REGISTER_ALL_COTES_INTEGRATE_METHODS_FOR_TYPE(FloatTests, float, float,
        float_configs[0]
    );

    // Можно также регистрировать отдельные методы для специфичных тестов
    REGISTER_ALL_GAUS_INTEGRATE_METHODS_FOR_TYPE(GAUS_N_POINT, double, double,
        variant_configs[0],
        variant_configs[1],
        variant_configs[2],
        variant_configs[3],
        variant_configs[4],
        variant_configs[5]
    );
}

namespace Integral_Computing_Function {

}

// ============================================================================
// ==================== ТЕСТЫ ДЛЯ Matrix2D ====================================
// ============================================================================
#include "../cpp_rewriting_SPBGU/Matrix2D.cpp"

// Подключаем функцию multiply (она в том же файле, что и Matrix2D)
// Функция multiply уже определена в Matrix2D.cpp

// ============================================================================
// Все тесты Matrix2D в одном пространстве имён
// ============================================================================
namespace Matrix2Dtest {

// ============================================================================
// Группа 1: Тесты конструкторов и деструктора
// ============================================================================
namespace Constructors {

    TEST(Matrix2DConstructorsTest, DefaultConstructor) {
        Matrix2D<double> m;
        EXPECT_EQ(m.getcol(), 0);
        EXPECT_EQ(m.getrow(), 0);
        EXPECT_EQ(m.getStorageOrder(), Matrix2D<double>::StorageOrder::RowMajor);
    }

    TEST(Matrix2DConstructorsTest, ConstructorWithSize) {
        Matrix2D<double> m(3, 4);
        EXPECT_EQ(m.getcol(), 3);
        EXPECT_EQ(m.getrow(), 4);
        EXPECT_EQ(m.getStorageOrder(), Matrix2D<double>::StorageOrder::RowMajor);
    }

    TEST(Matrix2DConstructorsTest, ConstructorWithSizeAndStorage) {
        Matrix2D<double> m(3, 4, Matrix2D<double>::StorageOrder::ColumnMajor);
        EXPECT_EQ(m.getcol(), 3);
        EXPECT_EQ(m.getrow(), 4);
        EXPECT_EQ(m.getStorageOrder(), Matrix2D<double>::StorageOrder::ColumnMajor);
    }

    TEST(Matrix2DConstructorsTest, SquareMatrixConstructor) {
        Matrix2D<double> m(5);
        EXPECT_EQ(m.getcol(), 5);
        EXPECT_EQ(m.getrow(), 5);
    }

    TEST(Matrix2DConstructorsTest, InitializerListConstructor) {
        Matrix2D<double> m = {{1, 2, 3}, {4, 5, 6}};
        EXPECT_EQ(m.getcol(), 3);
        EXPECT_EQ(m.getrow(), 2);
        EXPECT_DOUBLE_EQ(m(0, 0), 1);
        EXPECT_DOUBLE_EQ(m(0, 1), 2);
        EXPECT_DOUBLE_EQ(m(0, 2), 3);
        EXPECT_DOUBLE_EQ(m(1, 0), 4);
        EXPECT_DOUBLE_EQ(m(1, 1), 5);
        EXPECT_DOUBLE_EQ(m(1, 2), 6);
    }

    TEST(Matrix2DConstructorsTest, CopyConstructor) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2(m1);
        EXPECT_EQ(m2.getcol(), 2);
        EXPECT_EQ(m2.getrow(), 2);
        EXPECT_DOUBLE_EQ(m2(0, 0), 1);
        EXPECT_DOUBLE_EQ(m2(1, 1), 4);
    }

    TEST(Matrix2DConstructorsTest, MoveConstructor) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2(std::move(m1));
        EXPECT_EQ(m2.getcol(), 2);
        EXPECT_EQ(m2.getrow(), 2);
        EXPECT_DOUBLE_EQ(m2(0, 0), 1);
    }

} // namespace Constructors

// ============================================================================
// Группа 2: Тесты операторов доступа и геттеров/сеттеров
// ============================================================================
namespace Accessors {

    TEST(Matrix2DAccessorsTest, GetColGetRow) {
        Matrix2D<double> m(5, 7);
        EXPECT_EQ(m.getcol(), 5);
        EXPECT_EQ(m.getrow(), 7);
    }

    TEST(Matrix2DAccessorsTest, GetStorageOrder) {
        Matrix2D<double> m(2, 2, Matrix2D<double>::StorageOrder::ColumnMajor);
        EXPECT_EQ(m.getStorageOrder(), Matrix2D<double>::StorageOrder::ColumnMajor);
    }

} // namespace Accessors

// ============================================================================
// Группа 3: Тесты операторов класса (арифметика, присваивание, сравнение)
// ============================================================================
namespace Operators {

    // --- Операторы присваивания ---
    TEST(Matrix2DOperatorsTest, CopyAssignmentOperator) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2;
        m2 = m1;
        EXPECT_EQ(m2.getcol(), 2);
        EXPECT_EQ(m2.getrow(), 2);
        EXPECT_DOUBLE_EQ(m2(0, 0), 1);
        EXPECT_DOUBLE_EQ(m2(1, 1), 4);
    }

    TEST(Matrix2DOperatorsTest, MoveAssignmentOperator) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2;
        m2 = std::move(m1);
        EXPECT_EQ(m2.getcol(), 2);
        EXPECT_EQ(m2.getrow(), 2);
        EXPECT_DOUBLE_EQ(m2(0, 0), 1);
    }

    // --- Унарные операторы ---
    TEST(Matrix2DOperatorsTest, UnaryMinus) {
        Matrix2D<double> m = {{1, -2}, {3, -4}};
        Matrix2D<double> result = -m;
        EXPECT_DOUBLE_EQ(result(0, 0), -1);
        EXPECT_DOUBLE_EQ(result(0, 1), 2);
        EXPECT_DOUBLE_EQ(result(1, 0), -3);
        EXPECT_DOUBLE_EQ(result(1, 1), 4);
    }

    // --- Бинарные арифметические операторы ---
    TEST(Matrix2DOperatorsTest, Addition) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{5, 6}, {7, 8}};
        Matrix2D<double> result = m1 + m2;
        EXPECT_DOUBLE_EQ(result(0, 0), 6);
        EXPECT_DOUBLE_EQ(result(0, 1), 8);
        EXPECT_DOUBLE_EQ(result(1, 0), 10);
        EXPECT_DOUBLE_EQ(result(1, 1), 12);
    }

    TEST(Matrix2DOperatorsTest, Subtraction) {
        Matrix2D<double> m1 = {{5, 6}, {7, 8}};
        Matrix2D<double> m2 = {{1, 2}, {3, 4}};
        Matrix2D<double> result = m1 - m2;
        EXPECT_DOUBLE_EQ(result(0, 0), 4);
        EXPECT_DOUBLE_EQ(result(0, 1), 4);
        EXPECT_DOUBLE_EQ(result(1, 0), 4);
        EXPECT_DOUBLE_EQ(result(1, 1), 4);
    }

    TEST(Matrix2DOperatorsTest, Multiplication) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{5, 6}, {7, 8}};
        Matrix2D<double> result = m1 * m2;
        EXPECT_DOUBLE_EQ(result(0, 0), 19);
        EXPECT_DOUBLE_EQ(result(0, 1), 22);
        EXPECT_DOUBLE_EQ(result(1, 0), 43);
        EXPECT_DOUBLE_EQ(result(1, 1), 50);
    }

    TEST(Matrix2DOperatorsTest, MultiplicationRectangular) {
        // Матрица 4×3
        Matrix2D<double> m1 = {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9},
            {10, 11, 12}
        };
        // Матрица 3×4
        Matrix2D<double> m2 = {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12}
        };
        
        Matrix2D<double> result = m1 * m2;
        
        // Результат должен быть 4×4
        EXPECT_EQ(result.getrow(), 4);
        EXPECT_EQ(result.getcol(), 4);
        
        // Проверяем несколько элементов
        // result(0, 0) = 1*1 + 2*5 + 3*9 = 1 + 10 + 27 = 38
        EXPECT_DOUBLE_EQ(result(0, 0), 38);
        // result(0, 1) = 1*2 + 2*6 + 3*10 = 2 + 12 + 30 = 44
        EXPECT_DOUBLE_EQ(result(0, 1), 44);
        // result(0, 2) = 1*3 + 2*7 + 3*11 = 3 + 14 + 33 = 50
        EXPECT_DOUBLE_EQ(result(0, 2), 50);
        // result(0, 3) = 1*4 + 2*8 + 3*12 = 4 + 16 + 36 = 56
        EXPECT_DOUBLE_EQ(result(0, 3), 56);
        
        // result(1, 0) = 4*1 + 5*5 + 6*9 = 4 + 25 + 54 = 83
        EXPECT_DOUBLE_EQ(result(1, 0), 83);
        // result(1, 1) = 4*2 + 5*6 + 6*10 = 8 + 30 + 60 = 98
        EXPECT_DOUBLE_EQ(result(1, 1), 98);
        
        // result(3, 3) = 10*4 + 11*8 + 12*12 = 40 + 88 + 144 = 272
        EXPECT_DOUBLE_EQ(result(3, 3), 272);
    }

    TEST(Matrix2DOperatorsTest, MultiplyFunctionTest) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{5, 6}, {7, 8}};
        Matrix2D<double> result(2, 2);
        
        multiply(&m1, &m2, &result);
        
        EXPECT_DOUBLE_EQ(result(0, 0), 19);
        EXPECT_DOUBLE_EQ(result(0, 1), 22);
        EXPECT_DOUBLE_EQ(result(1, 0), 43);
        EXPECT_DOUBLE_EQ(result(1, 1), 50);
    }

    TEST(Matrix2DOperatorsTest, ScalarMultiplication) {
        Matrix2D<double> m = {{1, 2}, {3, 4}};
        Matrix2D<double> result = m * 2.0;
        EXPECT_DOUBLE_EQ(result(0, 0), 2);
        EXPECT_DOUBLE_EQ(result(0, 1), 4);
        EXPECT_DOUBLE_EQ(result(1, 0), 6);
        EXPECT_DOUBLE_EQ(result(1, 1), 8);
    }

    TEST(Matrix2DOperatorsTest, ScalarMultiplicationLeft) {
        Matrix2D<double> m = {{1, 2}, {3, 4}};
        Matrix2D<double> result = 2.0 * m;
        EXPECT_DOUBLE_EQ(result(0, 0), 2);
        EXPECT_DOUBLE_EQ(result(0, 1), 4);
    }

    TEST(Matrix2DOperatorsTest, ScalarDivision) {
        Matrix2D<double> m = {{2, 4}, {6, 8}};
        Matrix2D<double> result = m / 2.0;
        EXPECT_DOUBLE_EQ(result(0, 0), 1);
        EXPECT_DOUBLE_EQ(result(0, 1), 2);
        EXPECT_DOUBLE_EQ(result(1, 0), 3);
        EXPECT_DOUBLE_EQ(result(1, 1), 4);
    }

    // --- Операторы сравнения ---
    TEST(Matrix2DOperatorsTest, EqualityOperator) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{1, 2}, {3, 4}};
        EXPECT_TRUE(m1 == m2);
    }

    TEST(Matrix2DOperatorsTest, InequalityOperatorDifferentSize) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{1, 2, 3}, {4, 5, 6}};
        EXPECT_FALSE(m1 == m2);
    }

    TEST(Matrix2DOperatorsTest, InequalityOperatorDifferentValues) {
        Matrix2D<double> m1 = {{1, 2}, {3, 4}};
        Matrix2D<double> m2 = {{1, 2}, {3, 5}};
        EXPECT_FALSE(m1 == m2);
        EXPECT_TRUE(m1 != m2);
    }

    TEST(Matrix2DOperatorsTest, MultiplyAVX_FunctionTest) {
        // Создаём матрицы 4×4 для проверки правильности
        // A в RowMajor, B в ColumnMajor
        Matrix2D<double> A(4, 4, Matrix2D<double>::StorageOrder::RowMajor);
        Matrix2D<double> B(4, 4, Matrix2D<double>::StorageOrder::ColumnMajor);
        
        // Заполняем A
        A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3; A(0, 3) = 4;
        A(1, 0) = 5; A(1, 1) = 6; A(1, 2) = 7; A(1, 3) = 8;
        A(2, 0) = 9; A(2, 1) = 10; A(2, 2) = 11; A(2, 3) = 12;
        A(3, 0) = 13; A(3, 1) = 14; A(3, 2) = 15; A(3, 3) = 16;
        
        // Заполняем B
        B(0, 0) = 1; B(0, 1) = 2; B(0, 2) = 3; B(0, 3) = 4;
        B(1, 0) = 5; B(1, 1) = 6; B(1, 2) = 7; B(1, 3) = 8;
        B(2, 0) = 9; B(2, 1) = 10; B(2, 2) = 11; B(2, 3) = 12;
        B(3, 0) = 13; B(3, 1) = 14; B(3, 2) = 15; B(3, 3) = 16;
        
        Matrix2D<double> result(4, 4);
        
        multiplyAVX(&A, &B, &result);
        
        // Ожидаемые значения (скалярное умножение)
        // result(0, 0) = 1*1 + 2*5 + 3*9 + 4*13 = 1 + 10 + 27 + 52 = 90
        EXPECT_DOUBLE_EQ(result(0, 0), 90);
        // result(0, 1) = 1*2 + 2*6 + 3*10 + 4*14 = 2 + 12 + 30 + 56 = 100
        EXPECT_DOUBLE_EQ(result(0, 1), 100);
        // result(1, 0) = 5*1 + 6*5 + 7*9 + 8*13 = 5 + 30 + 63 + 104 = 202
        EXPECT_DOUBLE_EQ(result(1, 0), 202);
        // result(3, 3) = 13*4 + 14*8 + 15*12 + 16*16 = 52 + 112 + 180 + 256 = 600
        EXPECT_DOUBLE_EQ(result(3, 3), 600);
    }

    TEST(Matrix2DOperatorsTest, MultiplyAVX_PerformanceTest) {
        // Создаём большие матрицы для замера скорости
        const int M = 800;  // строки A
        const int K = 1200;  // столбцы A = строки B
        const int N = 800;  // столбцы B

        // A в RowMajor, B в ColumnMajor
        Matrix2D<double> A = Matrix2D<double>::random(K, M, 0.0, 1.0, Matrix2D<double>::StorageOrder::RowMajor);
        Matrix2D<double> B = Matrix2D<double>::random(N, K, 0.0, 1.0, Matrix2D<double>::StorageOrder::ColumnMajor);
        Matrix2D<double> result(N, M);

        // Замер времени AVX версии
        auto start = std::chrono::high_resolution_clock::now();
        multiplyAVX(&A, &B, &result);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double duration_ms = duration.count() / 1000.0;

        // Вывод результатов через cout (виден в логах Visual Studio)
        std::cout << std::endl;
        std::cout << "=== AVX Performance Results ===" << std::endl;
        std::cout << "Matrix Size: " << M << "x" << K << "x" << N << std::endl;
        std::cout << "A Storage: RowMajor" << std::endl;
        std::cout << "B Storage: ColumnMajor" << std::endl;
        std::cout << "AVX Time: " << duration.count() << " us (" << duration_ms << " ms)" << std::endl;
        std::cout << "===============================" << std::endl;

        // Сохраняем результаты через RecordProperty (для XML отчёта)
        RecordProperty("AVX_Time_us", duration.count());
        RecordProperty("AVX_Time_ms", duration_ms);
        RecordProperty("Matrix_Size", std::to_string(M) + "x" + std::to_string(K) + "x" + std::to_string(N));

        // Проверяем, что результат правильного размера
        EXPECT_EQ(result.getrow(), M);
        EXPECT_EQ(result.getcol(), N);

        // Проверяем, что результат не пустой
        bool nonZero = false;
        for (uint64_t i = 0; i < M && !nonZero; ++i) {
            for (uint64_t j = 0; j < N && !nonZero; ++j) {
                if (std::abs(result(i, j)) > 1e-10) {
                    nonZero = true;
                }
            }
        }
        EXPECT_TRUE(nonZero) << "Result matrix contains only zeros!";
    }

    TEST(Matrix2DOperatorsTest, MultiplyScalar_PerformanceTest) {
        // Создаём большие матрицы для замера скорости (те же размеры что и для AVX)
        const int M = 800;  // строки A
        const int K = 1200;  // столбцы A = строки B
        const int N = 800;  // столбцы B

        // A в RowMajor, B в ColumnMajor (для честного сравнения с AVX)
        Matrix2D<double> A = Matrix2D<double>::random(K, M, 0.0, 1.0, Matrix2D<double>::StorageOrder::RowMajor);
        Matrix2D<double> B = Matrix2D<double>::random(N, K, 0.0, 1.0, Matrix2D<double>::StorageOrder::ColumnMajor);
        Matrix2D<double> result(N, M);

        // Замер времени скалярной версии
        auto start = std::chrono::high_resolution_clock::now();
        multiply(&A, &B, &result);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double duration_ms = duration.count() / 1000.0;

        // Вывод результатов через cout (виден в логах Visual Studio)
        std::cout << std::endl;
        std::cout << "=== Scalar Performance Results ===" << std::endl;
        std::cout << "Matrix Size: " << M << "x" << K << "x" << N << std::endl;
        std::cout << "A Storage: RowMajor" << std::endl;
        std::cout << "B Storage: ColumnMajor" << std::endl;
        std::cout << "Scalar Time: " << duration.count() << " us (" << duration_ms << " ms)" << std::endl;
        std::cout << "=================================" << std::endl;

        // Сохраняем результаты через RecordProperty (для XML отчёта)
        RecordProperty("Scalar_Time_us", duration.count());
        RecordProperty("Scalar_Time_ms", duration_ms);
        RecordProperty("Matrix_Size", std::to_string(M) + "x" + std::to_string(K) + "x" + std::to_string(N));

        // Проверяем, что результат правильного размера
        EXPECT_EQ(result.getrow(), M);
        EXPECT_EQ(result.getcol(), N);

        // Проверяем, что результат не пустой
        bool nonZero = false;
        for (uint64_t i = 0; i < M && !nonZero; ++i) {
            for (uint64_t j = 0; j < N && !nonZero; ++j) {
                if (std::abs(result(i, j)) > 1e-10) {
                    nonZero = true;
                }
            }
        }
        EXPECT_TRUE(nonZero) << "Result matrix contains only zeros!";
    }

} // namespace Operators

} // namespace Matrix2Dtest

int main(int argc, char** argv) {
#if 1
    ::testing::GTEST_FLAG(catch_exceptions) = false;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
#else
    
#endif
}