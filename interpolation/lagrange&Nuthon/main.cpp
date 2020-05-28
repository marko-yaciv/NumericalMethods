#include <iostream>
#include <vector>
//#include <math.h>
using namespace std;
void fillXY(vector<double>& x, vector<double>& y);
double L(double xc,const vector<double>& x, const vector<double>& y);
double N(double xc,const vector<double>& x, const vector<double>& y);
double findDif(int n, vector<double> x, vector<double> y);
int main(){
    cout.precision(4);
    vector<double> x(10,0);
    vector<double> y(10,0);
    double xc = 1.461;
    fillXY(x,y);
    cout<<"Lagranja method"<<endl;
    double r = L(xc,x,y);
    cout<<"f("<<xc<<") = "<<r<<endl<<endl;

    cout<<"Nuthon Method"<<endl;
    r = N(xc,x,y);
    cout<<"f("<<xc<<") = "<<r<<endl;
    return 0;
}
void fillXY(vector<double>& x, vector<double>& y){
    x[0] = 1.430; y[0] = 0.880;
    x[1] = 1.435; y[1] = 0.889;
    x[2] = 1.440; y[2] = 0.890;
    x[3] = 1.445; y[3] = 0.891;
    x[4] = 1.450; y[4] = 0.892;
    x[5] = 1.455; y[5] = 0.893;
    x[6] = 1.460; y[6] = 0.894;
    x[7] = 1.465; y[7] = 0.895;
    x[8] = 1.470; y[8] = 0.896;
    x[9] = 1.475; y[9] = 0.897;

}
double findDif(int n ,vector<double> x, vector<double> y){
   double s = 0;
    for (int i = 0; i < n; ++i) {
        double d = 1;
        for (int j = 0; j < n; ++j) {
            if(j != i){
                d *= (x[i] - x[j]);
            }
        }
        double st = y[i]/d;
        s += st;
    }
    return s;
}
double L(double xc,const vector<double>& x, const vector<double>& y){
    double Ch;
    double Zn;
    double L=0;
    for (int i = 0; i < x.size(); i++) {
        Ch = 1; Zn = 1;
        for (int k = 0; k < x.size(); k++ ) {
            if ( k == i ) continue;
            Ch *= xc - x[k];
        }
        for(int k = 0; k < y.size();k++) {
            if (x[i] == x[k]) continue;
            Zn *= x[i] - x[k];
        }
        L += y[i] * Ch / Zn;
        cout << "\tL[" <<i+1 << "] = " << L << endl;
    }

    return L;

}
double N(double xc,const vector<double>& x, const vector<double>& y){
    vector<double> dif;
    dif.reserve(10);
    dif[0] = y[0];
    for (int i = 1; i < 10; ++i) {
        dif.push_back(findDif(i, x, y));
        cout<<"\tf(x"<<0<<"...x"<<i<<") = "<<dif.back()<<endl;
    }
    double s = 0;
    for (int i = 0; i < 10; ++i) {
        if(i == 0){
            s+= dif[0];
            continue;
        }
        double fuct = 1;
        for (int j = 0; j < i; ++j) {
            fuct *=(xc - x[j]);
        }
        fuct *= dif[i];
        s+=fuct;
    }
    return s;
}