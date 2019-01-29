//
//  main.cpp
//  IE531 MP2
//
//  Created by apple on 2/25/18.
//  Copyright © 2018 UIUC. All rights reserved.
//

//
//  main.cpp
//  Johnson_Lindenstrauss
//
//  Created by Ramavarapu Sreenivas on 4/3/17.
//  Copyright Â© 2017 Ramavarapu Sreenivas. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include <chrono>
#include <string.h>
#include <sstream>
#include "/Users/apple/Desktop/IE531/newmat11/newmatap.h"
#include "/Users/apple/Desktop/IE531/newmat11/newmat.h"
#include "/Users/apple/Desktop/IE531/newmat11//newmatio.h"


typedef std::pair<double,double> Box_Muller_Pair;

using namespace std;

int reduced_dimension, original_dimension;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

// U.I.I.D. RV generator
double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

// Using the Gaussian generator that is part of the C++ STL
double get_gaussian(double mean, double standard_deviation)
{
    std::normal_distribution<double> distribution(mean, standard_deviation);
    double number = distribution(generator);
    return (number);
}

int main(int argc, const char * argv[])
{
    int original_dimension, reduced_dimension, no_of_cols, no_of_trials;
    double epsilon, delta, diff;
    clock_t time_before, time_after;
    
    sscanf (argv[1], "%d", &original_dimension);
    sscanf (argv[2], "%d", &reduced_dimension);
    sscanf (argv[3], "%d", &no_of_cols);
    sscanf (argv[4], "%lf", &epsilon);
    sscanf (argv[5], "%lf", &delta);
    ifstream input_file(argv[6]);
    sscanf (argv[7], "%d", &no_of_trials);
    
    cout << "Johnson-Lindenstrauss Lemma Demo" << endl;
    cout << "Reading a (" << original_dimension << " x " << no_of_cols << ") Matrix from file '";
    cout << argv[6] << "'" << endl;
    cout << "Reduced Dimension = " << reduced_dimension << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "delta = " << delta << endl;
    cout << "Reduced Dimension (i.e. " << reduced_dimension << ") should be >= " << ceil(1/(epsilon*epsilon)*log(((double) no_of_cols)/delta));
    cout << " for the bound to hold with probability " << 1-delta << endl;
    
    Matrix A(original_dimension, no_of_cols);
    
    time_before = clock(); // recording time before we started to read data

    string x;
    string row;
    int i = 1;
    int j = 1;
    
    while(getline(input_file,row)){
        stringstream ss(row);
        while(getline(ss,x,',')){
        if(j <= original_dimension)
        {
            if (i <= no_of_cols) {
                A(j,i) = stod(x);
                i++;}
            if(i > no_of_cols )
            {
                j++;
                i = 1;
            }
        }
    }
    }

    time_after = clock(); // recording time after testing is complete
    diff = ((double) time_after - (double) time_before);
    
    cout << "It took " << (diff/CLOCKS_PER_SEC)/60.0 << " minutes to read data from file '" << argv[6] << "'" << endl;
    

    Matrix Generate_Random_Projection(reduced_dimension, original_dimension);

    for (int j = 1; j <= reduced_dimension; j++)
    {
        for (int i = 1; i <= original_dimension; i++)
        {
            Generate_Random_Projection(j,i) = get_gaussian(0, 1);
           // cout<< Generate_Random_Projection(j,i);
        }
    }

    
    RowVector N = Generate_Random_Projection.row(1);
    
    for (int j = 1; j <= reduced_dimension; j++)
    {
        for (int i = 1; i <= original_dimension; i++)
        {
        Generate_Random_Projection(j,i) = (sqrt(original_dimension/reduced_dimension))*Generate_Random_Projection(j,i)/NormFrobenius(N);
        }
    }

    
    
    Matrix R = Generate_Random_Projection;
    
    //cout<< R<<endl;

    // testing Johnson-Lindenstrauss Lemma
    int no_of_hits = 0;
    cout << "#Trails for the testing-phase = " << no_of_trials << endl;
    // this is the reduced-dimension representation of the x's (i.e the matrix of y's)
    
    
    Matrix C(reduced_dimension, no_of_cols);
    C = R*A;
    
   //cout<< C<<endl;
    
    time_before = clock(); // recording time before the testing starts
    
    // write code here for verification of JL-Lemma
    
    for(int i=1;i <= no_of_trials; i++)
        
    {
        int s1 = ceil(no_of_cols*get_uniform());
        int s2 = ceil(no_of_cols*get_uniform());
       
       ColumnVector E1 =A.column(s1);
       ColumnVector E2= A.column(s2);
        
       ColumnVector E3 = E2 -E1;
        
       ColumnVector F1 = R*E1;
       ColumnVector F2 = R*E2;
        
       ColumnVector F3 = F2 - F1;
        
        double lower = (1-epsilon)*NormFrobenius(E3);
        double upper =(1+epsilon)*NormFrobenius(E3);
        
       
      //cout<< F1 <<endl;
        
        if ((lower <= NormFrobenius(F3)) && (upper >= NormFrobenius(F3)))
            {
             no_of_hits ++;
            }
    }
    

     //cout<<no_of_hits;
    
    time_after = clock(); // recording time after testing is complete
    
    
    diff = ((double) time_after - (double) time_before);
    
    
    cout << "It took " << (diff/CLOCKS_PER_SEC)/60.0 << " minutes for testing to be completed" << endl;
    cout << "Johnson-Lindenstrauss Lemma is satisfied " << no_of_hits << "-many times over ";
    cout << no_of_trials << " attempts" << endl;
    cout << "Empirical Probability = " << ((double) no_of_hits/no_of_trials) << endl;
    cout <<"Theory says it should be at least "<< 1-delta << endl;

    return 0;
}

