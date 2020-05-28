#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

long double half_segment_method(long double LEnd, long double REnd, long double eps, int * cnt);
long double chord_method(long double LEnd,long double REnd,long double eps, int *cnt);
long double calcFuncValue(long double x);
long double calcChordValue(long double x);
long double calcDeriv(long double x);

int main() {
    long double LEnd, REnd,eps;
    ifstream InFile("source.txt");
    ofstream OutFile("result.txt");
    InFile >> LEnd >> REnd >> eps;
    int cnt = 0;
    long double DichotMethod = half_segment_method(LEnd, REnd, eps, &cnt);
    OutFile << fixed << "Dichotomia method result: " << DichotMethod << endl;
    OutFile <<"Iterations: "<< cnt <<"\n\n\n";

    cnt = 0;
    InFile >> LEnd >> REnd >> eps;
    long double ChordRes = chord_method(LEnd,REnd,eps, &cnt);
    if(ChordRes != -1){
        OutFile <<"Chord method result: "<< ChordRes <<endl;
        OutFile <<"iterations: " << cnt << endl;
    }else{
        OutFile << "There are no roots in this segment" << endl;
    }
    return 0;
}

long double half_segment_method(long double LEnd, long double REnd, long double eps, int * cnt){
    long double function_value;
    long double half_of_segment = 0;
    cout<<"L "<<LEnd<<"\t R " << REnd<<" ";

    while(abs(LEnd - REnd) >= eps){
        half_of_segment = (LEnd + REnd)/2.0;
        function_value = calcFuncValue(LEnd) * calcFuncValue(half_of_segment);

        if(function_value < 0){
            REnd = half_of_segment;
        }else{
            LEnd = half_of_segment;
        }
        *cnt += 1;
    }
    return half_of_segment;
}

long double chord_method(long double LEnd,long double REnd,long double eps, int *cnt){
    long double arg_value = -1;
    long double epsylon;

    epsylon = abs(LEnd - REnd);
    if(calcChordValue(LEnd) * calcChordValue(REnd) < 0){
        if(calcChordValue(REnd) * calcDeriv(REnd) > 0){
            while(epsylon >= eps){
                arg_value = LEnd - ( (calcChordValue(LEnd)*(REnd-LEnd)) /
                            (calcChordValue(REnd) - calcChordValue(LEnd)) );
                epsylon = abs(arg_value - LEnd);
                LEnd  = arg_value;
                *cnt +=1 ;
            }
        }else if(calcChordValue(LEnd) * calcDeriv(LEnd) > 0){
            while(epsylon >= eps){
                arg_value = REnd - ( (calcChordValue(REnd)*(REnd-LEnd)) /
                            (calcChordValue(REnd) - calcChordValue(LEnd)) );
                epsylon = abs(arg_value - LEnd);
                LEnd  = arg_value;
                *cnt +=1 ;
            }
        }
    }
    return arg_value;
}

long double calcFuncValue(long double x){
    return  1.8*pow(x, 2) - sin(10.0 * x);
}
long double calcChordValue(long double x){
    return (x * x * x) + (3 * x * x) + (6 * x) - 1;
}
long double calcDeriv(long double x){
    return (3*x*x) + (6*x) + 6;
}

