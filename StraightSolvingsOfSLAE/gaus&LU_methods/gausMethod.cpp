#include <iostream>
#include <fstream>
#define ROWS 3
#define COLS 3
using namespace std;
void fillMatr(double*** pKoef);
void fillResVect(double**pRes);
double findDet(double** matr);
double * gauss(double **pKoef, double *pRes);

double* gaus2(double **pKoef, double *pRes){
    double ext[ROWS][COLS+1];
    double max = 0;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS+1; ++j) {
            if(j < COLS){
                ext[i][j] = pKoef[i][j];
            }
            if(j == COLS){
                ext[i][j] = pRes[i];
            }
        }
    }

    max = ext[0][0];
    int ind = 0;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS+1; ++j) {
            if(ext[i][j] > max) {
                max = ext[i][j];
                ind = i;
            }
        }
    }
    cout<<"MAX: "<<max<<"   IND: "<<ind<<endl<<"M:"<<endl;
    double m[ROWS+1] = {0};
    for (int i = 0; i < ROWS+1; ++i) {
        for (int k = 0; k < ROWS + 1; ++k) {
            if(i!=ind){
                m[k] = ext[i][k]/ext[ind][k];
                ext[i][k] = ext[i][k] - ext[ind][k] * m[k];
            }
        }
    }
    for (double i : m) {
        cout<<i<<" ";
    }
    cout<<endl<<"---------------"<<endl;
    for (int i = 0; i < ROWS+1; ++i) {
        if(i != ind){
            for (int j = 0; j < ROWS + 1; ++j) {

            }
        }
    }
    for (int l = 0; l < ROWS + 1; ++l) {
        for (int i = 0; i < ROWS + 1; ++i) {
            cout<<ext[l][i]<<"\t";
        }
        cout<<endl;
    }

}

ofstream out("GausMethod.txt");
int main() {
    auto** pKoef = new double*[ROWS];
    for(int i=0;i<ROWS;++i){
        pKoef[i] = new double[COLS];
    }
    auto* pRes = new double[ROWS];
    auto * pAnswer = new double[ROWS];

    fillMatr(&pKoef);
    fillResVect(&pRes);
    cout << "Matrix: " << endl;
    for (int i = 0; i < ROWS; ++i) {
        for (int j = 0; j < COLS; ++j) {
            cout << pKoef[i][j]<<"\t";
        }
        cout << " = " << pRes[i];
        cout << endl;
    }
    double oDet = findDet(pKoef);
    if(!oDet){
        exit(1);
    }
    pAnswer = gauss(pKoef,pRes);
    for(int i=0;i<ROWS;++i)
        cout<<"x["<<i<<"]="<<pAnswer[i]<<endl;
    return 0;
}

void fillMatr(double*** pKoef){
//    (*pKoef)[0][0] = 0.43;
//    (*pKoef)[0][1] = 1.24;
//    (*pKoef)[0][2] = -0.58;
//    (*pKoef)[1][0] = 0.74;
//    (*pKoef)[1][1] = 0.83;
//    (*pKoef)[1][2] = 1.17;
//    (*pKoef)[2][0] = 1.43;
//    (*pKoef)[2][1] = -1.58;
//    (*pKoef)[2][2] = 0.83;
    (*pKoef)[0][0] =10;
    (*pKoef)[0][1] = 1;
    (*pKoef)[0][2] = -1;
    (*pKoef)[1][0] = 1;
    (*pKoef)[1][1] =10;
    (*pKoef)[1][2] =-1;
    (*pKoef)[2][0] =-1;
    (*pKoef)[2][1] =1;
    (*pKoef)[2][2] =10;
}
void fillResVect(double**pRes){
//    (*pRes)[0] = 2.71;
//    (*pRes)[1] = 1.26;
//    (*pRes)[2] = 1.03;
    (*pRes)[0] =11;
    (*pRes)[1] =10;
    (*pRes)[2] =10;
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

double * gauss(double **pKoef, double *pRes){
    cout.precision(4);
    double *x, max;
    int k, index;
    const double eps = 0.00001;
    x = new double[ROWS];
    k = 0;
    while (k < ROWS)
    {
        // Finding the row with the highest pKoef[i][k]
        max = abs(pKoef[k][k]);
        index = k;
        for (int i = k + 1; i < ROWS; i++)
        {
            if (abs(pKoef[i][k]) > max)
            {
                max = abs(pKoef[i][k]);
                index = i;
            }
        }
        //Swapping rows
        if (max < eps)
        {
            cout << "Can't get the solve because of non-zero column ";
            cout << index << " of A matrix" << endl;
            return 0;
        }
        for (int j = 0; j < ROWS; j++)
        {
            double temp = pKoef[k][j];
            pKoef[k][j] = pKoef[index][j];
            pKoef[index][j] = temp;
        }
        double temp = pRes[k];
        pRes[k] = pRes[index];
        pRes[index] = temp;

        // Normalisation of equancies
        for (int i = k; i < ROWS; i++)
        {
            double temp = pKoef[i][k];
            if (abs(temp) < eps) continue; // for zero koeficient continue
            for (int j = 0; j < ROWS; j++)
                pKoef[i][j] = pKoef[i][j] / temp;
            pRes[i] = pRes[i] / temp;
            if (i == k)  continue;
            for (int j = 0; j < ROWS; j++)
                pKoef[i][j] = pKoef[i][j] - pKoef[k][j];
            pRes[i] = pRes[i] - pRes[k];
        }
        k++;
        cout<<" "<<endl;
        for (int i = 0; i < ROWS; ++i) {
            for (int j = 0; j < COLS; ++j) {
                cout << pKoef[i][j] << "\t";
            }
            cout <<"\t  "<<pRes[i]<< endl;
        }

    }
    for (k = ROWS - 1; k >= 0; k--)
    {
        x[k] = pRes[k];
        for (int i = 0; i < k; i++)
            pRes[i] = pRes[i] - pKoef[i][k] * x[k];
    }
    out<<endl;

    return x;
}