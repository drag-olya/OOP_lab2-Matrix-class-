// Created by Olha Drahomeretska, group k-23
//

#include <iostream>
#include <vector>

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

    void operator=(Matrix& a) {

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                elements[i * n + j] = a.elements[i * n + j];
            }
        }

    }

    bool operator==(Matrix& a) {

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (elements[i * n + j] != a.elements[i * n + j]) {
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
                b.elements[i * n + j] = elements[i * n + j] + x;
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
                b.elements[i * n + j] = elements[i * n + j] + a.elements[i * n + j];
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
                b.elements[i * n + j] = elements[i * n + j] - x;
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
                b.elements[i * n + j] = elements[i * n + j] - a.elements[i * n + j];
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
                b.elements[i * n + j] = elements[i * n +j] * x;
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
                    b.elements[i * n + j] += elements[n * i + k] * a.elements[j + k * m];
                }
            }
        }
        return b;
    }

    Matrix operator~() {

        Matrix b = set_size(b, n, m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
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
                a.elements[i * a.n + j] = elem;
            }
            
        }
        //return is;
    }

    friend void operator<<(ostream & os, Matrix  a) {
        for (int i = 0; i < a.n; ++i) {
            for (int j = 0; j < a.m; ++j) {
                //os << a.elements[i * a.n + a.m] << " ";
                if (a.elements[i * a.n + j] < 10) {
                    cout << a.elements[i * a.n + j] << "  ";
                }
                else {
                    cout << a.elements[i * a.n + j] << " ";
                }
                
            }
            cout << endl;
            //os << endl;
        }
        //return os;  /// ??
    }

    };




    int main()
    {

        ///////// 2 matrices for test
        Matrix test1 = test1.set_size(test1, 3, 3);
        Matrix test2 = test2.set_size(test2, 3, 3);
        for (int i = 0; i < test1.n; i++)
        {
            for (int j = 0; j < test1.m; j++)
            {
                test1.elements[i * test1.n + j] = i * test1.n + j;
            }
        }
        cout << test1;
        cout << endl;
    
        for (int i = 0; i < test2.n; i++)
        {
            for (int j = 0; j < test2.m; j++)
            {
                test2.elements[i * test2.n + j] = (i + j)*4;
            }
        }

        cout << test2;
        cout << endl;
        ///////////////////

        cout << test1 * test2 << endl;


        //Matrix a, b;

        //cin >> a;
        //cin >> b;

        //cout << a;
        //cout << b;
        //cout << endl;
        
        //Matrix c = a + 3;
        //cout << a + 3;

        return 0;
    }