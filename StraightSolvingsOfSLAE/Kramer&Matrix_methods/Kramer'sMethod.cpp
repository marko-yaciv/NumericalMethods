#include <iostream>
#include <fstream>
#define ROWS 3
#define COLS 3

void fillMatr(double*** pKoef);
void fillResVect(double**pRes);
double findDet(double** matr);
void KfindEachDet(double***pMatr, double **pRes, double Odet, double **pAnswer);
using namespace std;
ofstream out("KramerMetod.txt");
int main(){
    auto** pKoef = new double*[ROWS];
    for(int i=0;i<ROWS;++i){
        pKoef[i] = new double[COLS];
    }
    auto* pRes = new double[ROWS];
    auto * pAnswer = new double[ROWS];

    out<<"---Kramer's Method---"<<endl;

    fillMatr(&pKoef);
    fillResVect(&pRes);
    out << "Matrix: " << endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out << pKoef[i][j]<<"\t";
        }
        out << " = " << pRes[i];
        out << endl;
    }


    double oDet = findDet(pKoef);
    if(!oDet){
        exit(1);
    }
    out<<endl<<"Ordinary Dterminant = "<<oDet<<endl;
    KfindEachDet(&pKoef,&pRes,oDet,&pAnswer);

    out<<"Roots: "<<endl;
    for(int i=0;i<ROWS;++i) {
        out << "x" << i + 1 << " = " << pAnswer[i]<<endl;
    }

    for(int i=0;i<ROWS;++i){
        delete pKoef[i];
    }
    delete[] pKoef;
    delete[] pRes;
    delete[] pAnswer;
}
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
void fillResVect(double**pRes){
    (*pRes)[0] = 2.71;
    (*pRes)[1] = 1.26;
    (*pRes)[2] = 1.03;
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
void KfindEachDet(double***pMatr, double **pRes, double Odet, double **pAnswer){

    for(int i=0;i<3;++i){
        if(i==0){
            (*pMatr)[0][0] = (*pRes)[0];
            (*pMatr)[1][0] = (*pRes)[1];
            (*pMatr)[2][0] = (*pRes)[2];
            (*pAnswer)[0] = findDet(*pMatr);
            out<<"Det1 = "<<(*pAnswer)[0]<<endl;
        }
        if(i==1){
            (*pMatr)[0][1] = (*pRes)[0];
            (*pMatr)[1][1] = (*pRes)[1];
            (*pMatr)[2][1] = (*pRes)[2];
            (*pAnswer)[1]  = findDet(*pMatr);
            out<<"Det2 = "<<(*pAnswer)[1]<<endl;
        }
        if(i==2){
            (*pMatr)[0][2] = (*pRes)[0];
            (*pMatr)[1][2] = (*pRes)[1];
            (*pMatr)[2][2] = (*pRes)[2];
            (*pAnswer)[2]  = findDet(*pMatr);
            out<<"Det3 = "<<(*pAnswer)[2]<<endl;
        }
        fillMatr(pMatr);
    }
    (*pAnswer)[0] /= Odet;
    (*pAnswer)[1] /= Odet;
    (*pAnswer)[2] /= Odet;
}