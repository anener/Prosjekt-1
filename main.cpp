#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <ctime>

using namespace std;
void prosjektB(int, double, double*);
void prosjektC(int, double, double*, double*);
void ludcmp(double**, int, int*, double*);
void prosjektE(int, double, double*);

int main()
{
    int n;
    cout << "Enter an n value: ";
    cin >> n;
    double h = 1.0/(n+1.0);

    ofstream myfile;
    string k = "xu_n" + to_string(n) + ".txt";
    myfile.open(k);
    myfile << "n = " << n << "\n";
    myfile << "h = " << h << "\n";
    myfile << "i x[i] u[i]\n";

    double *x = new double[n];
    double *f = new double[n];
    double *u = new double[n];

    x[0] = 0.0;
    f[0] = 100.0*exp(0.0);
    u[0] = 0.0;
    myfile << 0 << " " << x[0] << " " << u[0] << "\n";
    for(int i=1; i < n; i++) {
        x[i] = h*i;
        f[i] = 100.0*exp(-10.0*x[i]);
        u[i] = 1.0 - (1-exp(-10.0))*x[i] - exp(-10.0*x[i]);
        myfile << i << " " << x[i] << " " << u[i] << "\n";
    }
    myfile.close();

    //prosjektB(n, h, f);
    //prosjektC(n, h, f, u);

    prosjektE(n, h, f);


    return 0;
}


void prosjektB(int n, double h, double *f) {
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *v = new double[n];
    double *w = new double[n];
    double *bt = new double[n];
    double *wt = new double[n];

    for(int i=0; i < n; i++) {
        a[i] = -1.0;
        b[i] = 2.0;
        c[i] = -1.0;
        w[i] = f[i]*(h*h);
    }
    bt[1] = b[1];
    wt[1] = w[1];

    clock_t start = clock();
    for(int i=2; i<n; i++) {
        bt[i] = b[i] - ((a[i-1]*c[i-1])/(bt[i-1]));
        wt[i] = w[i] - ((a[i-1]*wt[i-1])/(bt[i-1]));
    }

    v[n-1] = (wt[n-1])/(bt[n-1]);
    for(int i=n-2; i>0; --i) {
        v[i] = wt[i] - ((c[i]*v[i+1])/(bt[i]));
    }
    clock_t end = clock();
    double timeUsed = double(end - start)/CLOCKS_PER_SEC;
    ofstream myfile;
    string k = "prosjektB_n" + to_string(n) + ".txt";
    myfile.open(k);
    myfile << "n = " << n << "\n";
    myfile << "h = " << h << "\n";
    myfile << "tid = " << timeUsed << "\n";
    myfile << "i v[i]\n";
    for(int i=0; i<n; i++) {
        myfile << i << " "  << v[i] << "\n";
    }
    myfile.close();


    return;
}

void prosjektC(int n, double h, double *f, double* u) {
    double *b = new double[n];
    double *v = new double[n];
    double *w = new double[n];
    double *bt = new double[n];
    double *wt = new double[n];

    for(int i=0; i < n; i++) {
        b[i] = 2.0;
        w[i] = f[i]*(h*h);
    }
    bt[1] = b[1];
    wt[1] = w[1];

    clock_t start = clock();
    for(int i=2; i<n; i++) {
        bt[i] = b[i] - ((1.0)/(bt[i-1]));
        wt[i] = w[i] - ((-wt[i-1])/(bt[i-1]));
    }

    v[n-1] = (wt[n-1])/(bt[n-1]);
    for(int i=n-2; i>0; --i) {
        v[i] = wt[i] - ((-v[i+1])/(bt[i]));
    }
    clock_t end = clock();
    double timeUsed = double(end - start)/CLOCKS_PER_SEC;

    double* e = new double[n];
    double y;
    for(int i=1; i<n; i++) {
        y = (v[i] - u[i])/u[i];
        e[i] = log10(abs(y));
    }

    ofstream myfile;
    string k = "prosjektC_n" + to_string(n) + ".txt";
    myfile.open(k);
    myfile << "n = " << n << "\n";
    myfile << "h = " << h << "\n";
    myfile << "tid = " << timeUsed << "\n";
    myfile << "relativ error = " << *max(e, e+n-1) << "\n";
    myfile << "i v[i]\n";
    for(int i=0; i<n; i++) {
        myfile << i << " "  << v[i] << "\n";
    }
    myfile.close();

    return;
}

void prosjektE(int n, double h, double *f) {
    double **a = new double*[n];
    for(int i=0; i<n; i++) {
        a[i] = new double[n];
    }
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            if(j==i) {
                a[i][j] = 2.0;
            }
            else if(j==(i+1) || j==(i-1)) {
                a[i][j] = -1.0;
            }
            else {
                a[i][j] = 0.0;
            }
        }
    }

    int *indx = new int[n];
    double *d = new double[n];
    double *y = new double[n];
    double *b = new double[n];
    double *v = new double[n];
    for(int i=0; i<n; i++) {
        b[i] = f[i]*(h*h);
    }


    ludcmp(a, n, indx, d);

    clock_t start = clock();
    y[1] = b[1]/a[1][1];
    for(int i=2; i<n; i++) {
        y[i] = b[i] - (a[i][i-1]*y[i-1]);
    }

    v[n-1] = (y[n-1])/(a[n-1][n-1]);
    for(int i=n-2; i>0; --i) {
        v[i] = (1/a[i][i]) * (y[i] - (a[i][i+1]*v[i+1]));
    }
    clock_t end = clock();
    double timeUsed = double(end - start)/CLOCKS_PER_SEC;

    ofstream myfile;
    string k = "prosjektE_n" + to_string(n) + ".txt";
    myfile.open(k);
    myfile << "n = " << n << "\n";
    myfile << "h = " << h << "\n";
    myfile << "tid = " << timeUsed << "\n";
    myfile << "i v[i]\n";
    for(int i=0; i<n; i++) {
        myfile << i << " "  << v[i] << "\n";
    }
    myfile.close();

}

void ludcmp(double **a, int n, int *indx, double *d)
{
   int      i, imax, j, k;
   double   big, dum, sum, temp, *vv;

  vv = new(nothrow) double [n];
  if(!vv) {
    printf("\n\nError in function ludcm():");
    printf("\nNot enough memory for vv[%d]\n",n);
    exit(1);
  }

   *d = 1.0;                              // no row interchange yet
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = 0.0;
      for(j = 0; j < n; j++) {
         if((temp = fabs(a[i][j])) > big) big = temp;
      }
      if(big == 0.0) {
         printf("\n\nSingular matrix in routine ludcmp()\n");
         exit(1);
      }
      vv[i] = 1.0/big;                 // save scaling */
   } // end i-loop */

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i< j; i++) {   // not i = j
         sum = a[i][j];
     for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
     a[i][j] = sum;
      }
      big = 0.0;   // initialization for search for largest pivot element
      for(i = j; i< n; i++) {
         sum = a[i][j];
     for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
     a[i][j] = sum;
     if((dum = vv[i]*fabs(sum)) >= big) {
        big = dum;
        imax = i;
     }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0;k< n; k++) {       // yes
        dum        = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k]    = dum;
     }
     (*d)    *= -1;            // and change the parit of d
     vv[imax] = vv[j];         // also interchange scaling factor
      }
      indx[j] = imax;
      if(fabs(a[j][j]) < 0.0)  a[j][j] = 0.0;

        /*
        ** if the pivot element is zero the matrix is singular
        ** (at least to the precision of the algorithm). For
        ** some application of singular matrices, it is desirable
        ** to substitute ZERO for zero,
        */

      if(j < (n - 1)) {                   // divide by pivot element
         dum = 1.0/a[j][j];
     for(i=j+1;i < n; i++) a[i][j] *= dum;
      }
   } // end j-loop over columns

   delete [] vv;   // release local memory

}
