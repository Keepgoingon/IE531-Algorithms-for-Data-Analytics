//
//  main.cpp
//  IE531PM4
//
//  Created by apple on 4/10/18.
//  Copyright © 2018 UIUC. All rights reserved.
//

#include <iostream>

//
//  main.cpp
//  Multivariate Gaussian via Metropolis-Hastings
//
//  Created by Ramavarapu Sreenivas on 3/10/17.
//  Copyright Â© 2017 Ramavarapu Sreenivas. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>

#include "/Users/apple/Desktop/IE531/newmat11/newmatap.h"
#include "/Users/apple/Desktop/IE531/newmat11/newmat.h"
#include "/Users/apple/Desktop/IE531/newmat11//newmatio.h"

#define PI 3.141592654

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double get_gaussian(double mean, double standard_deviation)
{
    std::normal_distribution<double> distribution(mean, standard_deviation);
    double number = distribution(generator);
    return (number);
}

double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

ColumnVector Generate_Independent_Multivariate_Gaussian(ColumnVector mean)
{
    ColumnVector x(mean.nrows());
    // Write some C++ code here that will generate "x" that has a Multivariate Gaussian
    //  RV "x" that has a mean of "mean" and a matrix that is the identity matrix
    for (int j = 1; j <= mean.nrows(); j++)
    {
        x(j) = get_gaussian(mean(j),1);
    
    }

    return x;
}

double MH_Discriminant(ColumnVector Current_value, ColumnVector Previous_value, SymmetricMatrix C, ColumnVector mean)
{
    double result;
    RowVector diff1 = (Current_value - mean).t();
    RowVector diff2 = (Previous_value - mean).t();
    SymmetricMatrix C1 = C.i();
    Matrix h1 = -0.5*(diff1 * C1 * (Current_value - mean));
    //cout<<h1<<endl;
    double h11 = h1(1,1);
  
    Matrix h2 = -0.5*(diff2 * C1 * (Previous_value - mean));
    double h22 = h2(1,1);
    // cout<<h2<<endl<<endl;
    double d1 =exp(h11);
                   
    double d2 =exp(h22);
    
    result =d1/d2;
    //cout<<result<<endl<<endl;

    if (result >= 1 )
        return 1;
    else
    // Write some C++ code here that computes equation 2 of the assignment description
       //cout<<result<<endl;
    return (result);
}

double Theoretical_PDF(ColumnVector x, SymmetricMatrix C, ColumnVector mean)
{
    // write C++ code that computes the expression of equation 1 of the assignment description
    double det;
   // cout<<"C is "<<C<<endl;
    det = C.Determinant();
    //cout<<"det is "<<det<<endl;
    
    double num;
    num = 1/sqrt(pow(2*PI, mean.nrows())*det);
    RowVector diff3 =(x- mean).t();
    ColumnVector diff4 = x-mean;
    SymmetricMatrix C2 = C.i();
    Matrix h3= diff3 * C2 * diff4;
    //cout<<"matrix is "<<h3<<endl;
    double num2;
    num2 = exp(-0.5 * h3(1,1));
    double num3;
    num3 = num*num2;
   // cout<<"number is"<<h3(1,1)<<endl;
    return (num3);
}

int main (int argc, char* argv[])
{
    ColumnVector y_prev, y_current;
    Matrix count(100,100);
    int no_of_trials, dimension;
    
    // 2D case
    dimension = 2;
    
    sscanf (argv[1], "%d", &no_of_trials);
    ofstream pdf_data(argv[2]);
    ofstream pdf_theory(argv[3]);
    
    // The covariance matrix
    SymmetricMatrix C(2);
    C(1,1) = 1.0;
    C(1,2) = 0.5;
    C(2,1) = 0.5;
    C(2,2) = 1.0;
    
    // The mean vector
    ColumnVector mean(2);
    mean(1) = 1.0;
    mean(2) = 2.0;
    
    cout << "Multivariate Gaussian Generator using MCMC-MH" << endl;
    cout << "Dimension = " << mean.nrows() << endl;
    cout << endl << "Mean Vector = " << endl << mean;
    cout << endl << "Covariance Matrix = " << endl << C;
    
    for (int i = 1; i <= 100; i++)
        for (int j = 1; j <= 100; j++)
            count(i,j) = 0.0;
    
    y_prev = Generate_Independent_Multivariate_Gaussian(mean);
   
    for (int i = 0; i < no_of_trials; i++)
    {
        y_current = Generate_Independent_Multivariate_Gaussian(mean);
        //cout<<"y_current is "<< y_current<<endl;
        
   
        if (get_uniform() < MH_Discriminant(y_current, y_prev, C, mean))
        {
             //int num22;
            //num22=0;
             //cout<<"MH_Discriminant is "<<MH_Discriminant(y_current, y_prev, C, mean)<<endl;
           // cout<<"y_current(1) is" <<y_current(1)<<endl;
           // cout<<"y_current(2) is" <<y_current(2)<<endl;
            for (int j = 1; j <= 100; j++) {
                for (int k = 1; k <= 100; k++) {
                    if ( (y_current(1) >= ((double) (j-52)/10)) && (y_current(1) < ((double) (j-51)/10)) &&
                        (y_current(2) >= ((double) (k-52)/10)) && (y_current(2) < ((double) (k-51)/10)) )
                   
                        count(j,k)++;
                   //cout<<"count("<<j<<","<<k<<") is "<<count(j,k)<<endl;
                }
            }
            y_prev = y_current;
        }
    }
    
    for (int j = 1; j <= 100; j++) {
        for (int k = 1; k <= 100; k++) {
            if (k < 100)
                pdf_data << count(j,k)/((double) no_of_trials) << ", ";
            if (k == 100)
                pdf_data << count(j,k)/((double) no_of_trials) << endl;
        }
    }
    
    double x1, x2;
    for (int j = 1; j <= 100; j++) {
        x1 = ((double) (j-51)/10);
        for (int k = 1; k <= 100; k++) {
            x2 = ((double) (k-51)/10);
            ColumnVector x(2);
            x(1) = x1;
            x(2) = x2;
            if (k < 100)
                pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << ", ";
            if (k == 100)
                pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << endl;
        }
    }
}
