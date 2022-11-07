// Created by Olha Drahomeretska, group k-23
//

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Matrix {

public:

    int n, m;
    vector<float> elements;

    Matrix set_size(Matrix a, int n, int m) {

        a.n = n; a.m = m;
        a.elements.resize(a.n * a.m);
        return a;
    }

    void operator=(Matrix a) {

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                elements[i * m + j] = a.elements[i * m + j];
            }
        }

    }

    bool operator==(Matrix& a) {

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (elements[i * m + j] != a.elements[i * m + j]) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix operator+(float x) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                b.elements[i * m + j] = elements[i * m + j] + x;
            }
        }
        return b;
    } // ?????

    Matrix operator+(Matrix& a) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                b.elements[i * m + j] = elements[i * m + j] + a.elements[i * m + j];
            }
        }
        return b;
    }

    Matrix operator-(float x) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                b.elements[i * m + j] = elements[i * m + j] - x;
            }
        }
        return b;
    } // ?????

    Matrix operator-(Matrix& a) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                b.elements[i * m + j] = elements[i * m + j] - a.elements[i * m + j];
            }
        }
        return b;
    }

    Matrix operator*(float x) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                b.elements[i * m + j] = elements[i * m +j] * x;
            }
        }
        return b;
    }

    Matrix operator*(Matrix& a) {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    b.elements[i * m + j] += elements[n * i + k] * a.elements[j + k * m];
                }
            }
        }
        return b;
    }

    Matrix operator~() {

        Matrix b = set_size(b, m, n);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                b.elements[i * n + j] = elements[j * m + i];
            }
        }
        return b;
    }

    friend void operator>>(istream& is, Matrix& a) {

        cout << "Enter size of the matrix:" << endl;
        //is >> a.n;
        //is >> a.m;
        cin >> a.n;
        cin >> a.m;

        a.elements.resize(a.n * a.m);

        cout << "Enter elements of the matrix:" << endl;

        for (int i = 0; i < a.n; ++i) {
            for (int j = 0; j < a.m; ++j) {

                float elem;
                cin >> elem;
                a.elements[i * a.m + j] = elem;
            }
            
        }
        //return is;
    }

    friend void operator<<(ostream & os, Matrix  a) {

        for (int i = 0; i < a.n; i++) {
            for (int j = 0; j < a.m; j++) {
                //os << a.elements[i * a.n + a.m] << " ";
                if (a.elements[i * a.m + j] < 10) {
                    cout << a.elements[i * a.m + j] << "  ";
                }
                else {
                    cout << a.elements[i * a.m + j] << " ";
                }
                
            }
            cout << endl;
            //os << endl;
        }
        //return os;  /// ??
    }

    vector<float> operator[](int row) {

        vector<float> nrow;
        nrow.resize(m);

        for (int i = 0; i < m; i++)
        {
            nrow[i] = elements[row * m + i];
        }

        return nrow;
    }

    };

    Matrix REF(Matrix m);

    vector<float> calc_res(Matrix m);


// a. метод √ауса розвТ€занн€ —Ћј–

    vector<float> GaussianElimination(Matrix m) { // m - розширена матриц€ n*n+1

        vector<float> res;

        Matrix ref_mat = REF(m);

        if (ref_mat == m) {
            cout << "Matrix is singular";
        }
        else {
            res = calc_res(ref_mat);
        }

        return res;

    }

    void swap_rows(Matrix m, int i, int j) {

        for (int k = 0; k < m.m; k++)
        {
            float cur = m[i][k];
            m[i][k] = m[j][k];
            m[j][k] = cur;
        }
    }


// зведенн€ до р€дковоњ ступ≥нчастоњ форми (Row echelon form)
Matrix REF(Matrix m) {

    for (int i = 0; i < m.n; i++)
    {
        // знаходимо найб≥льший за модулем у стовпц≥ ≥

        int i_max = i;
        float v_max = m[i_max][i];

        for (int k = i + 1; k < m.n; k++)
        {
            if (abs(m[k][i]) > v_max) {
                v_max = m[k][i];
                i_max = k;
            }
        }

        if (not m[i][i_max]) {

            return m; // вироджена матриц€ (singular matrix)
        }

        // р€док з визначеним найб≥льшим значенн€м стаЇ поточним
        if (i_max != i) {
            swap_rows(m, i, i_max);
        }

        // в≥д нижн≥х р€дк≥в в≥дн≥маЇмо ведучий домножений на коеф≥ц≥Їнт x
        for (int k = i + 1; k < m.n; k++)
        {

            float x = m[k][i] / m[i][i];

            for (int j = i + 1; j < m.m; j++)
            {
                m.elements[k * m.m + j] = m[k][j] - m[i][j] * x;

            }

            // елементи п≥д головною д≥агоналлю заповнюЇмо нул€ми
            m.elements[k * m.m + i] = 0;
        }
        //cout << m; //  test
    }
    //cout << endl << m; // test

    return m;
}


vector<float> calc_res(Matrix m) {

    vector<float> res;
    res.resize(m.n);

    for (int i = m.n - 1; i > -1; i--)
    {
        res[i] = m[i][m.m - 1];

        for (int j = i + 1; j < m.n; j++)
        {
            res[i] -= res[j] * m[i][j];
        }

        res[i] /= m[i][i];
    }

    return res;
}


// b.ћетод якоб≥ дл€ знаходженн€ власних значень матриц≥ (метод обертанн€ якоб≥)

float max_j_offdiag_in_row_i(Matrix m, int i) {

    float j_max = i + 1;
    float v_max = m[i][j_max];

    for (int j = i + 2; j < m.n; j++) {

        if (m[i][j] > v_max) {
            j_max = j;
            v_max = m[i][j];
        }
    }
    return j_max;
}


vector<float> JacobiEigenvalueAlgorithm(Matrix A) { // m - симетрична матриц€ n*n

    for (int i = 0; i < A.n - 1; i++)
    {
        float j_max = max_j_offdiag_in_row_i(A, i);

        while (abs(A[i][j_max]) > 0.1) {

            j_max = max_j_offdiag_in_row_i(A, i);

            float teta;

            if (A[i][i] - A[j_max][j_max]) {

                teta = 0.5 * atan(2 * A[i][j_max] / (A[i][i] - A[j_max][j_max]));

            }
            else {
                teta = M_PI / 4;
            }

            Matrix J = J.set_size(J, A.n, A.m);

            // заповнюЇмо матриц€ обертанн€ (якоб≥)
            for (int J_i = 0; J_i < J.n; J_i++)
            {
                for (int J_j = 0; J_j < J.m; J_j++) {

                    // заповнюЇмо д≥агональ
                    if (J_i == J_j) {
                        if (J_i == i or J_i == j_max) {

                            J.elements[J_i * J.m + J_j] = cos(teta);
                        }
                        else {
                            J.elements[J_i * J.m + J_j] = 1;
                        }
                    }

                    else if (J_i == i and J_j == j_max) {
                        J.elements[J_i * J.m + J_j] = sin(teta);
                    }
                    else if (J_i == j_max and J_j == i) {
                        J.elements[J_i * J.m + J_j] = -sin(teta);
                    }
                    else {
                        J.elements[J_i * J.m + J_j] = 0;
                    }
                }
            }

            A = ~J * A * J;

        } 

        }

    cout << endl << A;

    vector<float> res;
    res.resize(A.n);
    for (int i = 0; i < A.n; i++)
    {
        res[i] = A[i][i];
    }

    return res;

    }


Matrix rand_sym_matr(int n) {

    Matrix test_M = test_M.set_size(test_M, n, n);

    srand(time(0));

    for (int i = 0; i < test_M.n; i++) {
        for (int j = i; j < test_M.m; j++) {
            test_M.elements[i * test_M.m + j] = test_M.elements[j * test_M.m + i] = (float)rand()/1000;
        }
    }

    return test_M;
}



    int main()
    {

        ///////// 2 matrices for test
        Matrix test1 = test1.set_size(test1, 3, 3);
        Matrix test2 = test2.set_size(test2, 3, 3);
        for (int i = 0; i < test1.n; i++)
        {
            for (int j = 0; j < test1.m; j++)
            {
                test1.elements[i * test1.m + j] = i * test1.m + j;
            }
        }
        cout << test1;
        cout << endl;
    
        for (int i = 0; i < test2.n; i++)
        {
            for (int j = 0; j < test2.m; j++)
            {
                test2.elements[i * test2.m + j] = (i + j)*4;
            }
        }

        cout << test2;
        cout << endl;

        ///////////////////

        //cout << test1 * test2;

        //test for GaussianElimination()
        Matrix g = g.set_size(g, 3, 4);
        g.elements = {3, 2, -4, 3, 2, 3, 3, 15, 5, -3, 1, 14};

        vector<float> test_gaus = GaussianElimination(g);
        for (int i = 0; i < g.n; i++)
        {
            cout << "x" << i << " = " << test_gaus[i] << endl;
        }
        //


        //test for JacobiEigenvalueAlgorithm()
        Matrix a = a.set_size(a, 3, 3);
        a.elements = { 7, 5, 2, 5, 1, 0, 2, 0, 3 };

        float sq2 = sqrt(2);
        //a.elements = { 1, sq2, 2, sq2, 3, sq2, 2, sq2, 1};

        Matrix test_Jac = rand_sym_matr(3);
        cout << test_Jac;

        vector<float> test_jac = JacobiEigenvalueAlgorithm(test_Jac);
        for (int i = 0; i < a.n; i++)
        {
            cout << "e" << i << " = " << test_jac[i] << endl;
        }

        //



        //Matrix a, b;

        //cin >> a;
        //cin >> b;

        //cout << a;
        //cout << endl << b;
        //cout << endl;


        //Matrix c = a + 3;
        //cout << a - b;
        //cout << ~(a - b);
        

        return 0;
    }