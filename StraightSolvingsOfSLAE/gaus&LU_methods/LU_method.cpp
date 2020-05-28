#include <iostream>
#include <fstream>
#define ROWS 3
#define COLS 3
using namespace std;
void fillMatr(double*** pKoef);
void fillResVect(double**pRes);
double findDet(double** matr);
double* ulMethod(double**pKoef, double*pRes);
ofstream out("LU_method.txt");
int main() {
    auto** pKoef = new double*[ROWS];
    for(int i=0;i<ROWS;++i){
        pKoef[i] = new double[COLS];
    }
    auto* pRes = new double[ROWS];
    auto * pAnswer = new double[ROWS];

    fillMatr(&pKoef);
    fillResVect(&pRes);
    out.precision(4);
    out << "Matrix: " << endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out << pKoef[i][j]<<"\t";
        }
        out << "  " << pRes[i];
        out << endl;
    }
    double oDet = findDet(pKoef);
    if(!oDet){
        exit(1);
    }

    pAnswer = ulMethod(pKoef,pRes);
    for(int i=0;i<ROWS;++i)
        cout<<"x["<<i<<"]="<<pAnswer[i]<<endl;
    return 0;
}
void fillMatr(double*** pKoef){
    (*pKoef)[0][0] = 3;
    (*pKoef)[0][1] = 3;
    (*pKoef)[0][2] = 0;
    (*pKoef)[1][0] = 0;
    (*pKoef)[1][1] = 2;
    (*pKoef)[1][2] = 2;
    (*pKoef)[2][0] = 1;
    (*pKoef)[2][1] = 0;
    (*pKoef)[2][2] = 1;
}
void fillResVect(double**pRes){
    (*pRes)[0] = -2;
    (*pRes)[1] = -1;
    (*pRes)[2] = 0;
}
double findDet(double** matr){
    double res;
    res = ((matr[0][0]*matr[1][1]*matr[2][2])+
           (matr[0][1]*matr[1][2]*matr[2][0])+
           (matr[1][0]*matr[2][1]*matr[0][2]))-
          ((matr[0][2]*matr[1][1]*matr[2][0])+
           (matr[0][1]*matr[1][0]*matr[2][2])+
           (matr[0][0]*matr[1][2]*matr[2][1]));
    return res;
}
double* ulMethod(double**pKoef, double*pRes){

    auto* ans = new double[ROWS];
    auto **L = new double*[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        L[i] = new double[COLS];
    }
    auto **U = new double*[ROWS];
    for (int i = 0; i < ROWS; ++i) {
            U[i] = new double[COLS];
    }
    auto *Y = new double[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        Y[i] = 0;
        ans[i] = 0;
    }
    cout.precision(4);
    cout<<"------pKoef------"<<endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            cout << fixed << pKoef[i][j] << "  ";
        }
        cout<<endl;
    }
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }
    for (int i = 0; i < ROWS; ++i) {
            U[i][i] = 1;
    }
    for (int i = 0; i < ROWS; ++i) {
        L[i][0] = pKoef[i][0];
        U[0][i] = pKoef[0][i] / L[0][0];
    }
    double sum ;
    for (int i = 1; i < ROWS; i++)
    {
        for (int j = 1; j < ROWS; j++)
        {
            if (i >= j) // L-matrix
            {
                sum = 0;
                for (int k = 0; k < j; k++)
                    sum += L[i][k] * U[k][j];

                L[i][j] = pKoef[i][j] - sum;
            }
            else // U-matrix
            {
                sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[i][k] * U[k][j];

                U[i][j] = (pKoef[i][j] - sum) / L[i][i];
            }
        }
    }

    cout<<"--L--"<<endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            cout<<L[i][j]<<"    ";
        }
        cout<<endl;
    }
    cout<<"--U--"<<endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            cout<<U[i][j]<<"   ";
        }
        cout<<endl;
    }
    Y[0] = pRes[0]/L[0][0];
    for (int i = 1; i < ROWS; ++i) {
        sum = 0;
        for (int k = 0; k < i; ++k) {
            sum+=L[i][k] * Y[k];
        }
        Y[i] = (pRes[i] - sum)/L[i][i];
    }
    for (int l = 0; l < 3; ++l) {
        cout<<"Y"<<l<<"="<<Y[l]<<"  ";
    }
    cout<<endl;
    ans[2] = Y[2]/U[2][2];
    ans[1] = (Y[1] - U[1][2] * ans[2]) / U[1][1];
    ans[0] = (Y[0] - U[0][1]*ans[1] - U[0][2]*ans[2]) / U[0][0];
    return ans;
}