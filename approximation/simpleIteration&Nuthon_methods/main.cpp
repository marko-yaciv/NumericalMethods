#include <iostream>
#include <cmath>
#include <ctime>
#include <windows.h>
#include <chrono>
#define ROWS 2
using namespace std;

inline double calcF1(double x, double y);
inline double calcF2(double x, double y);
inline double derivF1(double x, double y);
inline double derivF2(double x, double y);
inline double xFromF1(double y);
inline double yFromF2(double x);

double determ(double** Arr, int size);
void findMinor(double** pMatr, double*** minor);
void transpo_matr(double ***pMatr, int r, int c);
void allocMatr(double ***pMatr, int r, int c);
bool isConvergent(double x, double y);
void makeYacobiMatr(double*** pMatr, double x, double y);

void iterWay(double& x, double& y);
void NuthonWay(double& x, double& y);
int main() {
    double x = 0.7;
    double y = 0.6;

    if(!isConvergent(x,y)){
        cout<<"the function isn't convergent"<<endl;
        return 1;
    }
    cout<<"Iterations way:\n";
    auto start = chrono::high_resolution_clock::now();
    iterWay(x,y);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
    //cout <<"Elapsed: "<<duration.count() <<"  microosec"<< endl;
    cout<< "X = "<<x<<";\tY = "<<y<<endl;
    x = 0.7;
    y = 0.6;

    cout<<"\nNuthon way\n";
    start = chrono::high_resolution_clock::now();
    NuthonWay(x,y);
    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(stop - start);
    //cout <<"Elapsed: "<<duration.count() <<"  microosec"<< endl;
    cout<< "X = "<<x<<";\tY = "<<y<<endl;

    return 0;
}
inline double calcF1(double x, double y){
    return cos(y) + x - 1.5;
}
inline double calcF2(double x, double y){
    return 2.0*y - sin(x - 0.5) - 1.0;
}
inline double derivF1(double x, double y){
    return -sin(y) + 1.0;

}
inline double derivF2(double x, double y){
    return 2.0 - cos(x - 0.5);
}
inline double xFromF1(double y){
    return 1.5 - cos(y);
}
inline double yFromF2(double x){
    return (1 + sin(x - 0.5))/2;
}

void allocMatr(double ***pMatr, int r, int c){
    *pMatr = new double*[r];
    for (int i = 0; i < r; ++i) {
        (*pMatr)[i] = new double[c];
    }
}

bool isConvergent(double x, double y){
    bool ec1 = (0 +  sin(y)) <= x && x < 1;
    bool ec2 = (0.5 * cos(0.5 - x) + 0) <= y && y < 1;

    return ec1 && ec2;
}

double determ(double** Arr, int size)
{
    int i,j;
    double det=0;
    double** matr;
    if(size==1)
    {
        det=Arr[0][0];
    }
    else if(size==2)
    {
        det=Arr[0][0]*Arr[1][1]-Arr[0][1]*Arr[1][0];
    }
    else
    {
        matr=new double*[size-1];
        for(i=0;i<size;++i)
        {
            for(j=0;j<size-1;++j)
            {
                if(j<i)
                    matr[j]=Arr[j];
                else
                    matr[j]=Arr[j+1];
            }
            det+=pow((double)-1, (i+j))*determ(matr, size-1)*Arr[i][size-1];
        }
        delete[] matr;
    }
    return det;
}

void transpo_matr(double ***pMatr, int r, int c) {
    double tmp;
    for (int i = 0; i < r; ++i) {
        for (int j = i; j < c; ++j) {
            tmp = (*pMatr)[i][j];
            (*pMatr)[i][j] = (*pMatr)[j][i];
            (*pMatr)[j][i] = tmp;
        }
    }
}
void findMinor(double** pMatr, double ***minor){
    (*minor)[0][0] = pMatr[1][1];
    (*minor)[0][1] = -pMatr[1][0];
    (*minor)[1][0] = -pMatr[0][1];
    (*minor)[1][1] = pMatr[0][0];
}
void makeYacobiMatr(double*** pMatr, double x, double y){
    (*pMatr)[0][0] = 1;
    (*pMatr)[0][1] = -sin(y);
    (*pMatr)[1][0] = -cos(-x + 0.5);
    (*pMatr)[1][1] = 2;
}

void iterWay(double& x, double& y){
    double eps = 10E-3;
    double xN , yN;

    while(true){
        xN = xFromF1(y);
        yN = yFromF2(x);
        cout<<"X = "<<x<<"   Y = "<<y<<endl;
        if(abs(xN - x) + abs(yN - y) < eps){
            break;
        }
        x = xN;
        y = yN;
    }
    x = xN;
    y = yN;

}
void NuthonWay(double& x, double& y){
    double eps = 10E-3;
    double **pMatr;
    allocMatr(&pMatr, 2, 2);
    cout.setf(ios_base::left);

    double xN;
    double yN;
    while(true){
        makeYacobiMatr(&pMatr, x,y);
        double det  = determ(pMatr, 2);
        cout<<"Det: "<<det<<endl;
        if(det!=0){
            double** minor;
            allocMatr(&minor, 2, 2);

            //wrapped  matrix <minor> to pMatr
            findMinor(pMatr, &minor);
            transpo_matr(&minor, 2, 2);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    minor[i][j] /= det;
                }
            }
            //delta of function in points x and y with negative sign
            auto* funcInPoint = new double[2];
            funcInPoint[0] = -calcF1(x,y);
            funcInPoint[1] = -calcF2(x,y);
            //finding next x and y;
            auto *delta = new double[2];
            delta[0] = (funcInPoint[0] * minor[0][0]) + (funcInPoint[1] * minor[1][0]);
            delta[1] = (funcInPoint[0] * minor[0][1]) + (funcInPoint[1] * minor[1][1]);

            cout << "delta X = " << delta[0] << "   delta Y = " << delta[1] << endl;
            xN = x + delta[0];
            yN = y + delta[1];
            if(abs(xN - x) + abs(yN- y) < eps){
                break;
            }
            x = xN;
            y = yN;
            cout<<"X = "<<x<<"   Y = "<<y<<endl;
        }else{
            cout<<"\nyacobi matr determinant  = 0\n";
            return;
        }
    }
    x = xN;
    y = yN;
    cout<<"X = "<<x<<"   Y = "<<y<<endl;
}

