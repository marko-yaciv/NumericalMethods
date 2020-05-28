#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> gauss(vector<vector<double> >& pKoef, vector<double>& pRes, const int ROWS){
    vector<double> x = {0,0,0,0};
    double max;
    int k, index;
    const double eps = 0.00001;
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
            exit(6);
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
            double tmp = pKoef[i][k];
            if (abs(tmp) < eps) continue; // for zero koeficient continue
            for (int j = 0; j < ROWS; j++)
                pKoef[i][j] = pKoef[i][j] / tmp;
            pRes[i] = pRes[i] / tmp;
            if (i == k)  continue;
            for (int j = 0; j < ROWS; j++)
                pKoef[i][j] = pKoef[i][j] - pKoef[k][j];
            pRes[i] = pRes[i] - pRes[k];
        }
        k++;
    }
    for (k = ROWS - 1; k >= 0; k--)
    {
        x[k] = pRes[k];
        for (int i = 0; i < k; i++)
            pRes[i] = pRes[i] - pKoef[i][k] * x[k];
    }

    return x;
}
vector<double> lowestSq(const vector<double>& x,const vector<double>& y, const int degree){
    const int s = degree+1;
    const int N = 5;
    int deg = 0;
    vector<vector<double> > pKoef(s,vector<double>(s));
    vector<double> pRes(s);
    vector<double> xd = {0,0,0,0,0,0,0};
    vector<double> a = {0,0,0,0};

    if(degree < 1 ||degree > 3){
        cout <<"Can solve only with 1 - 3 degree"<<endl;
        exit(0);
    }
    for (double& v:xd){
        for (int i = 0; i < N+1; i++)
            v += pow(x[i], deg);
        deg++;
    }
    cout<<"XD:\n";
    for (auto i : xd){
        cout<<i<<endl;
    }

//filling right side of the matrix
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            switch (i+j){
                    case 0:
                        pKoef[i][j] = N+1;
                        break;
                    case 1:
                        pKoef[i][j] = xd[1];
                        break;
                    case 2:
                        pKoef[i][j] = xd[2];
                        break;
                    case 3:
                        pKoef[i][j] = xd[3];
                        break;
                    case 4:
                        pKoef[i][j] = xd[4];
                        break;
                    case 5:
                        pKoef[i][j] = xd[5];
                        break;
                    case 6:
                        pKoef[i][j] = xd[6];
                        break;
                }
        }
    }
//filing left side of the matrix
    deg = 0;
    for (double & v : pRes) {
        for (int j = 0; j < N+1; ++j) {
            v += y[j] * pow(x[j], deg);
        }
        deg++;
    }
    cout.setf(ios_base::left);
//printing matrix
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            cout<<pKoef[i][j]<<"\t\t";
        }
        cout<<"| "<<pRes[i]<<endl;
    }
//solving system of equalities
    a = gauss(pKoef, pRes, s);
    return a;
}
int main() {
    vector<double> x = {8.03, 8.08, 8.16, 8.23, 8.26, 8.33};
    vector<double> y = {5.01, 4.78, 3.52, 3.12, 3.19, 2.95};

    vector<double> a = lowestSq(x,y,1);
    cout<<"P1(x) = "<<a[0]<<" + "<<a[1]<<"x"<<endl<<endl;
    a = lowestSq(x,y,2);
    cout<<"P2(x) = "<<a[0]<<" + "<<a[1]<<"x + "<<a[2]<<"x^2"<<endl<<endl;
    a = lowestSq(x,y,3);
    cout<<"P3(x) = "<<a[0]<<" + "<<a[1]<<"x + "<<a[2]<<"x^2 + "<<a[3]<<"x^3"<<endl;

    return 0;
}
