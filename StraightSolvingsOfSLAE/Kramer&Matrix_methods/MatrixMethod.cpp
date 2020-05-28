#include <iostream>
#include <cmath>
#include <fstream>
#define ROWS 3
#define COLS 3


double MfindDet(double** matr);
void MfindMinor(double**, double***);
void Mtranspo_matr(double ***pMatr);
void fillMatr(double***);

using namespace std;
int main() {

    auto **pKoef = new double *[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        pKoef[i] = new double[COLS];
    }
    auto *pRes = new double[ROWS];

    auto **pMinor = new double *[ROWS];
    for (int i = 0; i < ROWS; ++i) {
        pMinor[i] = new double[COLS];
    }

    pRes[0] = 2.71;
    pRes[1] = 1.26;
    pRes[2] = 1.03;
    ofstream out("MatrixMethodOut.txt");
    out << "---------Matrix Method---------" << endl<<endl;

    fillMatr(&pKoef);
    out << "   Matrix: " << endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out << pKoef[i][j]<<"\t";
        }
        out << " = " << pRes[i];
        out << endl;
    }

    double oDet = MfindDet(pKoef);
    if (!oDet) {
        exit(1);
    }
    out<<"Ordinaty determinant: "<<oDet<<endl<<endl;

    out << "---Minor matrix: ---" << endl;
    MfindMinor(pKoef, &pMinor);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out << pMinor[i][j]<<"\t";
        }
        out << endl;
    }

    out <<endl<< "----Transpoted minor matrix: ---" << endl;
    Mtranspo_matr(&pMinor);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out << pMinor[i][j]<<"\t";
        }
        out << endl;
    }
    //-----------wrapping matr---------------
    out <<endl<< "---Wrapped Matrix---" << endl;
    for (int i = 0; i < ROWS; ++i){
        for (int j = 0; j < COLS; ++j) {
            pKoef[i][j] = pMinor[i][j] / oDet;
            out << pKoef[i][j]<<"\t ";
        }
        out<<endl;
    }

    //-----finding result------------------------------
    auto * pAnswer = new double[ROWS];
    out<<"---Roots: ---"<<endl;
    for (int i = 0; i < ROWS; ++i){
        for (int j = 0; j < COLS; ++j){
            pKoef[i][j] *= pRes[j];
        }
        pAnswer[i] = pKoef[i][0]+ pKoef[i][1]+ pKoef[i][2];
        out << "x" << i + 1 << " = " << pAnswer[i]<<endl;
    }

    //----Free all----------------
    for(int i=0;i<ROWS;++i){
        delete[] pMinor[i];
        delete[] pKoef[i];
    }
    delete[] pMinor;
    delete[] pKoef;
    delete[] pRes;
    delete[] pAnswer;
}
//------------------------------------------------------------------------
void fillMatr(double*** pKoef){
    (*pKoef)[0][0] = 0.43;
    (*pKoef)[0][1] = 1.24;
    (*pKoef)[0][2] = -0.58;

    (*pKoef)[1][0] = 0.74;
    (*pKoef)[1][1] = 0.83;
    (*pKoef)[1][2] = 1.17;

    (*pKoef)[2][0] = 1.43;
    (*pKoef)[2][1] = -1.58;
    (*pKoef)[2][2] = 0.83;
}
//------------------------------------------------------------------------
double MfindDet(double** matr){
    double res;
    res = ((matr[0][0]*matr[1][1]*matr[2][2])+
           (matr[0][1]*matr[1][2]*matr[2][0])+
           (matr[1][0]*matr[2][1]*matr[0][2]))-
          ((matr[0][2]*matr[1][1]*matr[2][0])+
           (matr[0][1]*matr[1][0]*matr[2][2])+
           (matr[0][0]*matr[1][2]*matr[2][1]));
    return res;
}
//-------------------------------------------------------------------------
void MfindMinor(double **matr, double ***minor){
    (*minor)[0][0] = (matr[1][1]*matr[2][2])-(matr[2][1]*matr[1][2]);
    (*minor)[0][1] = -((matr[1][0]*matr[2][2])-(matr[2][0]*matr[1][2]));
    (*minor)[0][2] = (matr[1][0]*matr[2][1])-(matr[2][0]*matr[1][1]);

    (*minor)[1][0] = -((matr[0][1]*matr[2][2])-(matr[2][1]*matr[0][2]));
    (*minor)[1][1] = (matr[0][0]*matr[2][2])-(matr[2][0]*matr[0][2]);
    (*minor)[1][2] = -((matr[0][0]*matr[2][1])-(matr[2][0]*matr[0][1]));

    (*minor)[2][0] = (matr[0][1]*matr[1][2])-(matr[1][1]*matr[0][2]);
    (*minor)[2][1] = -((matr[0][0]*matr[1][2])-(matr[1][0]*matr[0][2]));
    (*minor)[2][2] = (matr[0][0]*matr[1][1])-(matr[1][0]*matr[0][1]);


}
//------------------------------------------------------------------------
void Mtranspo_matr(double ***pMatr) {
    double tmp;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = i; j < COLS; ++j) {
            tmp = (*pMatr)[i][j];
            (*pMatr)[i][j] = (*pMatr)[j][i];
            (*pMatr)[j][i] = tmp;
        }
    }
}
//-------------------------------------------------------------------------
