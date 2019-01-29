
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "/Users/apple/Desktop/IE531/newmat11/newmatap.h"
#include "/Users/apple/Desktop/IE531/newmat11/newmat.h"
#include "/Users/apple/Desktop/IE531/newmat11//newmatio.h"

using namespace std;

double Frobenious_Norm(Matrix Data)
{
    double ForbeniusNorm = 0;
    double ForbeniusNorm1 = 0;
   Data = Data * Data.t();
    for (int i = 1; i <= Data.nrows(); i++)
    {
        for (int j = 1; j <=Data.ncols(); j++)
        {
            ForbeniusNorm += Data(i,j) * Data(i,j);
        }
    }
    
    ForbeniusNorm1 = sqrt(ForbeniusNorm);
    return   ForbeniusNorm1;
}

double Norm(ColumnVector B1)
{
    double lnorm = 0;
    
    for (int i= 1; i<= B1.nrows(); i++)
    {
        lnorm += B1(i)*B1(i);
    }
    return lnorm;
}


Matrix Matrix_Sketch(Matrix Data, double epsilon)
{
    // Edo Liberty's Matrix Sketch will have the same number of rows
    // as the original data; the #cols is ceil(2.0/epsilon)
    int cols_of_sketch = ceil(2.0 / epsilon);
    
    if (cols_of_sketch < Data.nrows())
    {
        Matrix Result(Data.nrows(), cols_of_sketch);
        // write this part of the code yourself
        
        Matrix U, V;
        DiagonalMatrix D;
        Matrix minb;
        for (int i = 1; i <= Data.nrows(); i++)
        {
            for (int j = 1; j <= cols_of_sketch; j++)
            {
                Result(i, j) = Data(i, j);
            }
            
        }
        
        for (int j = cols_of_sketch + 1; j <= Data.ncols(); j++)
        {
            SVD(Result, D, U, V);
   
             Matrix minb = U*D*V;
            ColumnVector B1 =minb.column(cols_of_sketch);
            
            double delta = Norm(B1);
          for (int i = 1; i <= cols_of_sketch; i++)
            {
                if (D(i, i)*D(i, i) - delta> 0)
                {
                    D(i, i) = sqrt(D(i, i) * D(i, i) - delta);
                }
                else
                {
                    D(i, i) = 0;
                }
                
            }
            
            Result = U * D;
            
            for (int i = 1; i <= Data.nrows(); i++)
            {
                Result(i, cols_of_sketch) = Data(i, j);
            }
            
        }
        return Result;
    }
    else
    {
       Matrix Result(Data.nrows(),Data.nrows());
        
        Matrix U, V;
        DiagonalMatrix D;
        for (int i = 1; i <= Data.nrows(); i++)
        {
            for (int j = 1; j <= Data.nrows(); j++)
            {
                Result(i, j) = Data(i, j);
            }
            
        }
        
        for (int j = Data.nrows() + 1; j <= Data.ncols(); j++)
        {
            SVD(Result, D, U, V);
            
            Matrix minb = U*D*V;
            ColumnVector B1 =minb.column(Data.nrows());

             double delta = Norm(B1);

            for (int i = 1; i <= D.nrows(); i++)
            {
                if (D(i, i)*D(i, i) - delta> 0)
                {
                    D(i, i) = sqrt(D(i, i)*D(i, i) - delta);
                }
                else
                {
                    D(i, i) = 0;
                }
                
            }
            Result = U * D;
            for (int i = 1; i <= Data.nrows(); i++)
            {
                Result(i, Data.nrows()) = Data(i, j);
            }
        }
        return Result;
    }

    
}

int main(int argc, char* argv[])
{
    int dimension, no_of_data_points;
    double epsilon;
    
    sscanf(argv[1], "%d", &dimension);
    sscanf(argv[2], "%d", &no_of_data_points);
    sscanf(argv[3], "%lf", &epsilon);
    ifstream input_file(argv[4]);
    ofstream output_file(argv[5]);
    
    Matrix Data(dimension, no_of_data_points);
    
    cout << "Edo Liberty's Matrix Sketching Algorithm" << endl;
    cout << "----------------------------------------" << endl;
    cout << "Original Data-Matrix has " << dimension << "-rows & " << no_of_data_points << "-cols" << endl;
    cout << "Epsilon = " << epsilon << " (i.e. max. of " << 100 * epsilon << "% reduction of  Frobenius-Norm of the Sketch Matrix)" << endl;
    cout << "Input File = " << argv[4] << endl;
    
    // Read the Data
    for (int i = 1; i <= dimension; i++)
        for (int j = 1; j <= no_of_data_points; j++)
        {
            double x;
            input_file >> x;
            Data(i, j) = x;
        }
    
    // Compute the Frobenius-Norm of the original Data-Matrix
    double Data_Forbenius_Norm = Frobenious_Norm(Data);
    cout << "Frobenius Norm of the (" << Data.nrows() << " x " << Data.ncols() << ") Data Matrix = ";
    cout << Data_Forbenius_Norm << endl;
    
    Matrix Sketch(dimension, min(dimension, (int)ceil(2 / epsilon)));
    Sketch = Matrix_Sketch(Data, epsilon);
    double Sketch_Forbenius_Norm = Frobenious_Norm(Sketch);
    cout << "Frobenius Norm of the (" << Sketch.nrows() << " x " << Sketch.ncols() << ") Sketch Matrix = ";
    cout << Sketch_Forbenius_Norm << endl;
    cout << "Change in Frobenius-Norm between Sketch & Original  = ";
    cout << setprecision(3) << 100 * (Sketch_Forbenius_Norm - Data_Forbenius_Norm) / Data_Forbenius_Norm << "%" << endl;
    
    output_file << Sketch;
    cout << "File `" << argv[5] << "' contains a (" << Sketch.nrows() << " x " << Sketch.ncols();
    cout << ") Matrix-Sketch" << endl;
    
    
}
