#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
typedef long double ld;

inline ld IcalcFiValue(ld val);
inline ld IcalcFuncValue(ld val);
inline ld NcalcFuncValue(ld val);
inline ld NcalcDiverValue(ld val);
inline ld NcalcSecondDiverValue(ld val);
ld IcalcArgValue(ld LEnd, ld REnd, ld eps, int* cnt);
ld NcalcArgValue(ld LEnd, ld REnd, ld eps, int* cnt);
int main() {
    ld LEnd, REnd, eps;
    ifstream InFile("source.txt");
    ofstream OutFile("result.txt");
    InFile >> LEnd >> REnd >> eps;

    int cnt = 0;
    ld IterResult = IcalcArgValue(LEnd, REnd, eps, &cnt);
    OutFile<<"Iterations method: "<<endl;
    if(IterResult != -100){
        OutFile<<"Result: "<<IterResult<<endl;
        OutFile<<"Iterations: "<<cnt<<endl;
    }else{
        OutFile<<"There are NO roots in this range"<<endl;
    }
    cnt = 0;
    ld NutonResult = NcalcArgValue(LEnd,REnd,eps, &cnt);
    OutFile<<endl<<endl<<"Nuthon method:"<<endl;
    if(NutonResult!= -100){
        OutFile<<"Result: "<<NutonResult<<endl;
        OutFile<<"Iterations: "<<cnt<<endl;
    }else{
        OutFile<<"There are NO roots in this range"<<endl;
    }
    return 0;
}


inline ld IcalcFiValue(ld val){
    return 1/pow(2,val);
}
inline ld IcalcFuncValue(ld val){
    return val*pow(2,val) -1;
}
inline ld NcalcFuncValue(ld val){
    return (val*val*val) - 3*(val*val) + 9*val-10;
}
inline ld NcalcDiverValue(ld val){
    return 3*(val*val) - 6*val + 9;
}
inline ld NcalcSecondDiverValue(ld val){
    return 6*val - 6;
}

ld NcalcArgValue(ld LEnd, ld REnd, ld eps, int* cnt){
    ld curent, next, stop;
    if(NcalcFuncValue(LEnd)*NcalcFuncValue(REnd) < 0){
        if(NcalcFuncValue(LEnd)*NcalcSecondDiverValue(LEnd) > 0){
            curent = LEnd;
            do{
                next  = curent - (NcalcFuncValue(curent))/
                                 (NcalcDiverValue(curent));
                stop  = next - curent;
                curent = next;
                *cnt+=1;
            }while(abs(stop) >= eps);
        }else {
            curent = REnd;
            do {
                next = curent - (NcalcFuncValue(curent)) /
                                (NcalcDiverValue(curent));
                stop = next - curent;
                curent = next;
                *cnt+=1;
            }while (abs(stop) >= eps);
        }
        return curent;
    }else{
        return -100;
    }
}
ld IcalcArgValue(ld LEnd, ld REnd, ld eps, int *cnt){
    if(IcalcFuncValue(LEnd) * IcalcFuncValue(REnd) < 0){
        ld curent = LEnd;
        ld stop, next;
        do{
            next = IcalcFiValue(curent);
            stop = next-curent;
            curent = next;
            *cnt+=1;
        }while(abs(stop) >= eps);
        return curent;
    }else{
        return -100;
    }
}

