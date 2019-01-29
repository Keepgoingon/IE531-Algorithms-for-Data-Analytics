//
//  ie531pm42.cpp
//  IE531PM4
//
//  Created by apple on 4/10/18.
//  Copyright © 2018 UIUC. All rights reserved.
//

//
//  main.cpp
//  Multivariate Gaussian via Gibbs Sampling
//
//  Created by Ramavarapu Sreenivas on 3/8/17.
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

// this is just a hint -- if you are confused with what I suggest, you are welcome to take a
// different approach... no compulsion to use this thought-process, per se!



ColumnVector Gibbs_sampler_for_Multivariate_Gaussian(int index, ColumnVector Previous_value, SymmetricMatrix C, ColumnVector mean)
{
    if(index > mean.nrows())
     index = 1;
    else
        index = index;
    // we want the conditional distribution of the "index"-th rv, conditioned on the
    // values of the other indices (i.e. "index"-th rv is random; all others are fixed).
    // See desription and linked YouTube Video from Prof. Niemi.
    
    // Took the formulae from http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node7.html
    // Keep in mind that the expression for the correlation-matrix for the conditional probability
    // has a minor-typo at this reference (he has transposed "i" and "j"; you could have guessed it
    // by checking for dimensional consistency of the resulting correlation-matrix.
    
    // In the following -- I am assuming the original correlation matrix (Sigma) is partioned as
    // [Sigma11 Sigma12; Sigma21 Sigma22] (in MATLAB notation).  Sigma11 is an (n-1) x (n-1) matrix
    // Sigma12 is an (n-1)-Columnvector; Sigma21 is an (n-1)-Rowvector; Sigma22 is a scalar.
    // The import being -- we are keeping all the "top" (n-1)-many variables constant, and picking
    // a value for the bottom variable (i.e conditioning the bottom-variable on the previous (n-1)
    // variables.
    
    // Since the Gibbs-Sampler is not always going to sample the bottom-value, you need to go through
    // some clerical-work to make things work.
  
    // Write code to construct Sigma11 here
    Matrix Sigma11(mean.nrows()-1,mean.nrows()-1);
    for (int i=1; i <index; i++)
    {
        for (int j=1; j<index; j++)
        {
        Sigma11(i,j) = C(i,j);
        }
        for (int j=index+1; j<= mean.nrows(); j++)
        {
        Sigma11(i,j-1) = C(i,j);
        }
    }
    
       for (int i=index+1; i <= mean.nrows(); i++)
    {
        for (int j=1; j<index; j++)
        {
            Sigma11(i-1,j) = C(i,j);
        }
        for (int j=index+1; j<= mean.nrows(); j++)
        {
            Sigma11(i-1,j-1) = C(i,j);
        }
       
        Sigma11(index,index) = C(index,index);
    }
    
    //cout<<"Sigma11 is "<<Sigma11<<endl;
    // Write code to construct Sigma12 here
    
    ColumnVector Sigma12(mean.nrows()-1);
    for (int i=1; i < index; i++)
    {
        Sigma12(i) = C(i,index);
        
    }
    for (int i=index+1; i <= mean.nrows(); i++)
        {
        Sigma12(i-1) = C(i,index);
        }
   //  cout<<"Sigma12 is "<<Sigma12<<endl;
    // Write code to construct Sigma21 here
  
    RowVector Sigma21(mean.nrows()-1);
        for (int i=1; i < index; i++)
        {
            Sigma21(i) = C(index,i);
        }
        for (int i=index+1; i <= mean.nrows(); i++)
        {
            Sigma21(i-1) = C(index,i);
        }
       //  cout<<"Sigma21 is "<<Sigma21<<endl;
    
    // Write code to construct Sigma22 here
    
    Matrix Sigma22(1,1);
    
    Sigma22(1,1) = C(mean.nrows(),mean.nrows());
    
    //cout<<"Sigma22 is "<<Sigma22(1,1)<<endl;
        
    // Write some account-keeping code here... in my "teliology" https://www.merriam-webster.com/dictionary/teleology
    // xx is a ColumnVector, which starts off with xx = Previous_value; and then the "index"-th variable is an
    // appropriately generated Univariate-Normal RV (called "x")
    
    ColumnVector Present(mean.nrows());
    {
    Present = Previous_value;
    
       // cout <<"Present is "<< Present;
      
    double mean1;
    double cov;
      
       Matrix name[2][2];
       name[0][0] = Sigma11;
       name[0][1] = Sigma12;
       name[1][0] = Sigma21;
       name[1][1] = Sigma22;
        


        Matrix  mean11 = mean(index) + name[index-1][2-index] * name[2-index][2-index].i() * (Previous_value(3-index)-mean(3-index));
        
        mean1 = mean11(1,1);
       
     //   cout<< "mean1 is "<<mean1<<endl<<endl;
        
        Matrix cov11 = name[index-1][index-1] - name[2-index][index-1].t() * name[2-index][2-index].i() * name[2-index][index-1];
        cov = cov11(1,1);
        
       //  cout<< "cov is "<<mean1<<endl<<endl;
     Present(index) = get_gaussian(mean1, cov);
      //  cout<<"index is "<<index<<endl<<endl;
        
     //    cout<<"Present(index) is "<<Present(index)<<endl<<endl<<endl;
        
        index ++;
        
        
        //cout<<"present is "<<endl<<Present<<endl;
   
        return Present;
        }
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
    ColumnVector y_prev, y;
    Matrix count(100,100);
    int no_of_trials, dimension;
    
    // 2D case
    dimension = 2;
    
    sscanf (argv[1], "%d", &no_of_trials);
    ofstream pdf_data(argv[2]);
    ofstream pdf_theory(argv[3]);
    
    // The covariance matrix
    SymmetricMatrix C(2);
    C(1,1) = 0.75;
    C(1,2) = 0.25;
    C(2,1) = 0.25;
    C(2,2) = 0.5;
    
    // The mean vector
    ColumnVector mean(2);
    mean(1) = 1.0;
    mean(2) = 2.0;
    
    cout << "Multivariate Gaussian Generator using Gibbs Sampling" << endl;
    cout << "Dimension = " << mean.nrows() << endl;
    cout << endl << "Mean Vector = " << endl << mean;
    cout << endl << "Covariance Matrix = " << endl << C;
    
    for (int i = 1; i <= 100; i++)
        for (int j = 1; j <= 100; j++)
            count(i,j) = 0.0;
    
    y_prev = mean;
    for (int i = 0; i < no_of_trials; i++)
    {
        y = Gibbs_sampler_for_Multivariate_Gaussian(i%(mean.nrows())+1, y_prev, C, mean);
        for (int j = 1; j <= 100; j++) {
            for (int k = 1; k <= 100; k++) {
                if ( (y(1) >= ((double) (j-52)/10)) && (y(1) < ((float) (j-51)/10)) &&
                    (y(2) >= ((double) (k-52)/10)) && (y(2) < ((float) (k-51)/10)) )
                    count(j,k)++;
            }
        }
        y_prev = y;
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
