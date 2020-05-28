#include <iostream>
#include <fstream>
#include <cmath>
#define ROWS 5
#define COLS 3
using namespace std;
void fillMatr(double*** pKoef);
void fillResVect(double***);
void printMatr(double**pMatr, int, int);
void printMatr(double*pMatr, char);
void allocMatr(double***, int, int);
void allocMatr(double**, int);
void transpo_matr(double **pMatr, double **[], int r, int c);
void findSqRoots(double ***pMatr, double ***pRes, double ***pTrans, double **pAns);
ofstream out("lab6_out.txt");
int main(){
    double**pKoef;
    allocMatr(&pKoef,ROWS,COLS);
    fillMatr(&pKoef);

    double**pRes;
    allocMatr(&pRes,ROWS,1);
    fillResVect(&pRes);

    double **pTrans;
    allocMatr(&pTrans,COLS,ROWS);
    out<<"Matrix:\n";
    printMatr(pKoef, ROWS,COLS);
    out<<"\nTranspoted matrix:\n";
    transpo_matr(pKoef,&pTrans,ROWS,COLS);
    printMatr(pTrans,COLS,ROWS);

    double* pAns;
    allocMatr(&pAns,COLS);
    findSqRoots(&pKoef,&pRes,&pTrans,&pAns);
    out<<"-------X-------\n";
    printMatr(pAns,'X');

    for (int i = 0; i < ROWS; ++i) {
        delete[] pKoef[i];
    }
    for (int i = 0; i < COLS; ++i) {
        delete[] pTrans[i];
    }
    delete[] pKoef;
    delete[] pTrans;
    delete[] pRes;
    delete[] pAns;
    return 0;
}
void allocMatr(double ***pMatr, int r, int c){
    *pMatr = new double*[r];
    for (int i = 0; i < r; ++i) {
        (*pMatr)[i] = new double[c];
    }
}
void allocMatr(double** pVect, int s){
    *pVect = new double[s];
}
void fillResVect(double***pRes){
    (*pRes)[0][0] = 3;
    (*pRes)[1][0] = -3;
    (*pRes)[2][0] = -1;
    (*pRes)[3][0] = 6;
    (*pRes)[4][0] = 3;
}
void fillMatr(double*** pKoef){
    (*pKoef)[0][0] = 2;
    (*pKoef)[0][1] = -3;
    (*pKoef)[0][2] = 1;

    (*pKoef)[1][0] = -1;
    (*pKoef)[1][1] = -5;
    (*pKoef)[1][2] = 0;


    (*pKoef)[2][0] = 2;
    (*pKoef)[2][1] = -1;
    (*pKoef)[2][2] = 3;


    (*pKoef)[3][0] = 1;
    (*pKoef)[3][1] = -1;
    (*pKoef)[3][2] = -4;

    (*pKoef)[4][0] = -1;
    (*pKoef)[4][1] = 3;
    (*pKoef)[4][2] = 2;
}
void printMatr(double**pMatr,const int r,const int c){
    out.precision(3);
    out.setf(ios_base::right);

    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            out.width(3);
            out<<pMatr[i][j]<<"\t";
        }
        out<<endl;
    }
}
void printMatr(double*pMatr, char c){

    out.precision(3);
    for (int j = 0; j < COLS; ++j) {
        out<<c<<"["<<j<<"] = "<<pMatr[j]<<"\t";
    }
    out<<endl;
}
void transpo_matr(double **pMatr, double **trans[ROWS],const int r, const int c) {
    double tmp;
    for (int i = 0; i < c; ++i) {
        for (int j = 0; j < r; ++j) {
            (*trans)[i][j] = pMatr[j][i];
        }
    }

}
void findSqRoots(double ***pMatr, double ***pRes, double ***pTrans, double **pAns){
    double **N;
    allocMatr(&N,COLS,COLS);
    for (int i = 0; i < COLS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            for (int k = 0; k < ROWS; ++k) {
                N[i][j]  += (*pTrans)[i][k]*(*pMatr)[k][j];
            }
        }
    }
    out<<"N:\n";
    printMatr(N,COLS,COLS);

    double *C;
    allocMatr(&C,COLS);
    for (int i = 0; i < COLS; ++i) {
        for (int k = 0; k < ROWS; ++k) {
            C[i]  += (*pTrans)[i][k]*(*pRes)[k][0];
        }
    }
    out<<"C:\n";
    printMatr(C,'C');

    double **L;
    allocMatr(&L,COLS,COLS);
    L[0][1] = 0;
    L[0][0] = sqrt(N[0][0]);
    L[1][0] = N[1][0] / L[0][0];
    L[2][0] = N[2][0] / L[0][0];
    L[1][1] = sqrt(N[1][1] - pow(L[1][0],2));
    L[2][1] = (N[2][1] - L[2][0] * L[1][0]) / L[1][1];
    L[2][2] = sqrt(N[2][2] - pow(L[2][0],2) - pow(L[2][1],2));
    out<<"------L-----\n";
    printMatr(L,COLS,COLS);

    double **Lt;
    allocMatr(&Lt,COLS,COLS);
    transpo_matr(L,&Lt,COLS,COLS);
    out<<"------Lt-----\n";
    printMatr(Lt,COLS,COLS);

    double *Y;
    allocMatr(&Y,COLS);
    Y[0] = C[0]/L[0][0];
    Y[1] = (C[1] - L[1][0]*Y[0])/L[1][1];
    Y[2] = (C[2] - L[2][0]*Y[0] - L[2][1]*Y[1])/L[2][2];
    out<<"-------Y-------\n";
    printMatr(Y,'Y');

    (*pAns)[2] = Y[2] / Lt[2][2];
    (*pAns)[1] = (Y[1] - Lt[1][2] * (*pAns)[2])/Lt[1][1];
    (*pAns)[0] = (Y[0] - Lt[0][2] * (*pAns)[2] - Lt[0][1] * (*pAns)[1])/ Lt[0][0];

    for(int i = 0; i < COLS; ++i){
        delete[] N[i];
        delete[] L[i];
        delete[] Lt[i];
    }
    delete[] N;
    delete[] L;
    delete[] Lt;
    delete[] C;
    delete[] Y;
}
