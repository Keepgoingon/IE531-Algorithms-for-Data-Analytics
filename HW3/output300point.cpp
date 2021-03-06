
//
//  main.cpp
//  IE531PM3
//
//  Created by apple on 3/19/18.
//  Copyright © 2018 UIUC. All rights reserved.
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

using namespace std;
int exponent;
Matrix G ;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);
// U.I.I.D. RV generator


double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0, 1);
    double number = distribution(generator);
    return (number);
}

float Random()
{
    int a = -5;
    float num = (10*get_uniform() + a);
    return num;
}


Matrix repeated_squaring(Matrix B, int exponent )
{
    B(5,5);
    IdentityMatrix A(5);
    
    if (exponent == 0)
        return A;
    if (exponent == 1)
        return B;
    
    {
        if (exponent % 2 == 0)
            return(repeated_squaring(B*B, exponent/2));
        else
            return(B * repeated_squaring(B*B,(exponent-1)/2));
    }
}

Matrix direct_squaring(Matrix B, int exponent)
{
    Matrix G = B;
    for(int i=0; i<exponent-1;i++)
        G = B * G;
    return G;
}



int main(int argc, char * argv[])
{
    int exponent;
    double diff1,diff2;
    clock_t time_before1, time_after1,time_before2, time_after2;
    
    for(int i=0; i<=300; i++)
    {
        exponent= i;
    
    
    
    cout << exponent <<"      ";
    Matrix C(5,5);
    Matrix B(5,5);
    for (int j = 1; j <= 5; j++)
    {
        for (int i = 1; i <= 5; i++)
        {            B(j,i) = Random();
        }
    }
    //cout<<"Initial Matrix: "<<endl;
    
    
    // cout<<"Repeated Squaring Result: "<<endl;
    time_before1 = clock(); // recording time before we started to read data
    
    
    C = repeated_squaring (B, exponent);
    
    //cout << "Repeated Squaring Result:" << endl;
   // cout << setw(5)<<setprecision(2)<< C << endl;
    
    
    time_after1 = clock(); // recording time after testing is complete
    diff1 = ((double) time_after1 - (double) time_before1);
    
    cout << (diff1/CLOCKS_PER_SEC)/60.0 << "          ";
    
    //cout<<"Direct Squaring Result:"<<endl;
    
    time_before2 = clock();
    
    G = direct_squaring (B, exponent);
    
    // cout << "Direct Squaring Result : "<<endl;
    //cout <<setw(5)<<setprecision(2)<< G <<endl;
    
    time_after2 = clock();
    diff2 = ((float) time_after2 - (float) time_before2);
    cout << diff2/CLOCKS_PER_SEC <<endl;
    }
    
}










