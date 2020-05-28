#include <iostream>
#include <fstream>
#include <cmath>
#define ROWS 2
#define COLS 2
using namespace std;
void fillMatr(double*** pKoef);
void fillResVect(double**);
void printMatr(double**pMatr);
void printMatr(double*pMatr);
void allocMatr(double***, int, int);
void allocMatr(double**, int);
int yacobiMethod( double **pKoef,  double *pRes,double **ans);
int zeidelMethod( double **pKoef,  double *pRes,double **ans);
ofstream out("con");
int main() {
    double**pKoef;
    allocMatr(&pKoef,ROWS,COLS);
    fillMatr(&pKoef);

    double*pRes;
    allocMatr(&pRes,ROWS);
    fillResVect(&pRes);
    out<<"Matrix:"<<endl;
    printMatr(pKoef);
    out<<endl;
    int iterCnt = 0;
//--------------------------------------------------------
    auto *pAns = new double[ROWS];
    iterCnt = yacobiMethod(pKoef,pRes,&pAns);
    out.precision(3);
    out<<endl<<"<<<<<<<Yacobi's method>>>>>>>>>>>"<<endl;
    for (int l = 0; l < ROWS; ++l) {
        out << "X[" << l << "] = " << pAns[l] << "\n";
        pAns[l] = 0;
    }
    out<<"Iterations:"<<iterCnt<<endl;
    iterCnt = 0;
//--------------------------------------------------------
    out<<endl<<"<<<<<<Zeidel's method>>>>>>>>"<<endl;
    iterCnt = zeidelMethod(pKoef,pRes,&pAns);
    for (int l = 0; l < ROWS; ++l) {
        out << "X[" << l << "] = " << pAns[l] << "\n";
        pAns[l] = 0;
    }
    out<<"Iterations:"<<iterCnt<<endl;
    delete[] pKoef;
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

void fillResVect(double**pRes){
    (*pRes)[0] = -1;
    (*pRes)[1] = 15;
//    (*pRes)[2] = 2.34;
//    (*pRes)[3] = -0.72;
}
void fillMatr(double*** pKoef){
    (*pKoef)[0][0] = 5;
    (*pKoef)[0][1] = -3;
    //(*pKoef)[0][2] = -0.13;
    //(*pKoef)[0][3] = -0.11;

    (*pKoef)[1][0] = 1;
    (*pKoef)[1][1] = 7;
//    (*pKoef)[1][2] = 0.09;
//    (*pKoef)[1][3] = -0.06;
//
//    (*pKoef)[2][0] = 0.11;
//    (*pKoef)[2][1] = 0.05;
//    (*pKoef)[2][2] = -2.02;
//    (*pKoef)[2][3] = 0.12;
//
//    (*pKoef)[3][0] = 0.13;
//    (*pKoef)[3][1] = 0.1;
//    (*pKoef)[3][2] = 0.24;
//    (*pKoef)[3][3] = 0.43;
}

void printMatr(double**pMatr){
    out.setf(ios_base::showpos);
    out.setf(ios_base::right);
    out.precision(3);
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out<<pMatr[i][j]<<"\t";
        }
        out<<endl;
    }
}
void printMatr(double*pMatr){
    out.precision(3);
    for (int j = 0; j < COLS; ++j) {
        out<<pMatr[j]<<"   ";
    }
    out<<endl;
}

int yacobiMethod( double **pKoef, double *pRes, double **ans){
    const double eps = 0.001;
    cout.precision(3);
    out.setf(ios_base::showpos);
    out.setf(ios_base::left);
    double a[ROWS][COLS] = {0};
    double b[ROWS] = {0};
    //finding alpha and beta matrix
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            if(i!=j){
                a[i][j] = -pKoef[i][j]/pKoef[i][i];
            }
        }
        b[i] = pRes[i]/pKoef[i][i];
    }
    out<<"---a----:\n";
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            out<<a[i][j]<<"\t";
        }
        out<<endl;
    }
    out<<"----b----:\n";
    for (int j = 0; j < COLS; ++j) {
        out<<b[j]<<"\t";
    }out<<endl;
    //check for convergency
    for (int i = 0; i < ROWS; ++i) {
        int s = 0;
        for (int j = 0; j < COLS; ++j) {
            s+= abs(a[i][j]/a[i][i]);
        }
        if(s>1){
            out<<"Matrix isn't convergent"<<endl;
            exit(666);
        }
    }
    double prevAns[ROWS] = {0};
    for (int i = 0; i < ROWS; ++i) {
        prevAns[i] = (*ans)[i] = b[i];
    }
    double stop = 0;
    int cnt = 0;
    while(true){
        for (int i = 0; i < ROWS; ++i){
            double s1 = 0;
            for (int j = 0; j < i; ++j) {
                s1 += a[i][j] * prevAns[j];
            }
            double s2  = 0;
            for (int j = i+1; j < ROWS; ++j) {
                s2 += a[i][j]*prevAns[j];
            }
            (*ans)[i] = b[i] + s1 + s2;
            out<<"ans: "<<(*ans)[i]<<"  ";
        }
        out<<endl;
        double s3 = 0;
        for (int k = 0; k < ROWS; ++k) {
            s3 += pow((*ans)[k]-prevAns[k] ,2);
        }
        stop = sqrt(s3);
        if(stop < eps) break;
        for (int l = 0; l < ROWS; ++l) {
            prevAns[l] = (*ans)[l];
        }
        ++cnt;
    }
    return cnt;
}

int zeidelMethod( double **pKoef,  double *pRes,double **ans){
    const double eps = 0.001;
    cout.precision(3);
    double a[ROWS][COLS] = {0};
    double b[ROWS] = {0};
    //finding alpha and beta matrix
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            if(i!=j){
                a[i][j] = -pKoef[i][j]/pKoef[i][i];
            }
        }
        b[i] = pRes[i]/pKoef[i][i];
    }
    //check for convergency
    for (int i = 0; i < ROWS; ++i) {
        int s = 0;
        for (int j = 0; j < COLS; ++j) {
            s+= abs(a[i][j]/a[i][i]);
        }
        if(s>1){
            out<<"Matrix isn't convergent"<<endl;
            exit(666);
        }
    }
    double prevAns[ROWS] = {0};
    for (int i = 0; i < ROWS; ++i) {
        prevAns[i] = (*ans)[i] = b[i];
    }
    double stop = 0;int cnt = 0;
    while(true){
        for (int i = 0; i < ROWS; ++i){
            double s1 = 0;
            for (int j = 0; j < i; ++j) {
                s1 += a[i][j] * (*ans)[j];
            }
            double s2  = 0;
            for (int j = i+1; j < ROWS; ++j) {
                s2 += a[i][j]*prevAns[j];
            }
            (*ans)[i] = b[i] + s1 + s2;
        }
        double s3 = 0;
        for (int k = 0; k < ROWS; ++k) {
            s3 += pow((*ans)[k]-prevAns[k] ,2);
        }
        stop = sqrt(s3);
        if(stop < eps) break;
        for (int l = 0; l < ROWS; ++l) {
            prevAns[l] = (*ans)[l];
        }
        ++cnt;
    }
    return cnt;
}