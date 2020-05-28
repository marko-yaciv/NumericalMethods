#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

inline double funcInPoint(double x){
    return pow(x,3) * exp(x/3.0) * cosh(x);
}
double squareL(const double a, const double b, const double h){
    double vx = a;
    double res = 0;
    while(vx < b) {
        res += funcInPoint(vx);
        vx += h;
    }
    return res*h;
}
double squareR(const double a, const double b, const double h){
    double vx = a + h;
    double res = 0;
    while(vx <= b) {
        res += funcInPoint(vx);
        vx += h;
    }
    return res*h;
}
double squareC(const double a, const double b, const double h){
    double vx = a;
    double res = 0;
    while(vx < b) {
        res += funcInPoint(vx + h/2);
        vx += h;
    }
    return res*h;
}
double trapecia(const double a, const double b, const double h){
    double vx = a;
    double res = 0;
    while(vx < b) {
        res +=  (funcInPoint(vx) + funcInPoint(vx + h)) / 2;
        vx += h;
    }
    return res * h;
}
double simpson(const double a, const double b, const double h){
    double vx = a;
    double res = 0;
    while(vx <= b){
        res += funcInPoint(vx) + 4 * funcInPoint(vx + h) + funcInPoint(vx + 2 * h);
        vx += h * 2;
    }
    return (res * h)/3;
}
void printXY(const vector<double>& x,const vector<double>& y){
    auto sz = x.size();
    cout.precision(3);
    for(int i = 0; i < sz; ++i){
        cout <<i<<". "<< "x" <<" = " << x.at(i) <<"\t";
        cout << "y" << " = " << y.at(i) << endl;
    }
    cout<<endl;
}
int main() {
    double a = 0, b = 6;
    double eps = 10e-15;
    vector<double> x, y;

    cout << "Left rectangle method:" << endl;
    //calculating normal h
    double h = sqrt(eps);
    double h2 = h/2;
    while(abs(squareL(a,b,h) - squareL(a,b,h2)) > eps){
        h2 /= 2;
        h = h2;
    }
    cout << "I = " <<  squareL(a,b,h) << endl;

    cout << "Right rectangle method:" << endl;
    while(abs(squareR(a,b,h) - squareR(a,b,h2)) > eps){
        h2 /= 2;
        h = h2;
    }
    cout << "I = " << squareR(a,b,h) << endl;

    cout << "Center rectangle method:"<< endl;
    while(abs(squareC(a,b,h) - squareC(a,b,h2)) > eps){
        h2 /= 2;
        h = h2;
    }
    cout << "I = " << squareC(a,b,h) << endl;

    cout << "Trapecias method:" << endl;
    while(abs(trapecia(a,b,h) - trapecia(a,b,h2)) > eps){
        h2 /= 2;
        h = h2;
    }
    cout << "I = " << trapecia(a,b,h) << endl;


    cout << "Simpson method:" << endl;
    //calculating normal h
    h = pow(eps, 1.0/4.0);
    h2 = h/2;
    while(abs(simpson(a,b,h) - simpson(a,b,h)) > eps){
        h2 /= 2;
        h = h2;
    }
    cout << "I = " << simpson(a,b,h) << endl;
    return 0;
}
