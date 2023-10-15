#include <iostream>
#include "Header.h"
#include <fstream>
#include <string>
using namespace std;
int main()
{
    vector<double> a{ 0 , 0 };
    vector<double> b{ 1 , 0 };
    vector<double> c{ 0 , 1 };
    vector<double> d{ 1 , 1 };
    Way aw(a);
    Way bw(b);
    Way cw(c);
    Way dw(d);
    vector<double> q {0.0};
    vector<double> w {1.0};
    vector<double> e {1.0};
    vector<double> r {0.0};
    Way qw(q);
    Way ww(w);
    Way ew(e);
    Way rw(r);
    vector<Way> X  { aw, bw, cw, dw};
    vector<Way> Y { q, w, e, r };
    vector<int> razmer = { 2, 3, 1};
    Network setka(razmer);
    setka.ShowNetwork();
    setka.Train(X, Y, 0.5, 0.0000007, 10000, 0.7); // тренируем сеть
    setka.ShowNetwork();
    for (int i = 0; i < 4; i++)
    {
        Way output = setka.Forward(X[i]); //прогоняем вектор ренировочных значений
        cout << "X: ";
        cout << X[i].GetWay(0); //получаем первый эелемент
        cout << X[i].GetWay(1); // и второй
        cout << " Y: ";
        cout << Y[i].GetWay(0); // получаем правильное значение текущего тренировочного массива
        cout << " output: ";
        cout << output.GetWay(0);
        cout << "" <<endl;
    }
    string path = "Output_Matrix.txt";
    ofstream tempf;
    tempf.open(path);
    if (!tempf.is_open())
    {
        cout << "error of opening";
    }
    else
    {
        for (int k = 0; k < setka.layersN; k++)
        {
            for (int i = 0; i < razmer[k+1]; i++)
            {
                for (int j = 0; j < razmer[k]; j++)
                {
                    tempf << setka.weights[k].transpond(setka.weights[k]).GetMatrix(i, j)<< " ";
                    //tempf << setka.weights[k].GetMatrix(i, j);
                }
                tempf << endl;
            }
            tempf << " " << endl;
        }
        cout << "file has been opened" << endl;

    }
    return 0;
}
