#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <numeric>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>

using namespace std;

typedef vector<double> row;
typedef vector<row> matrix;
typedef vector<matrix> array3d;

void repmatRow(row * vec, int x, int y, int z, array3d *arr);
void repmatArr(array3d * arr, int x, int y, int z, array3d *result);
void AAminus(array3d * arr1, array3d * arr2, array3d *result);
void AAdiv(array3d * arr1, array3d * arr2, array3d *result);
void AAprod(array3d * arr1, array3d * arr2, array3d *result);
void AAsum(array3d * arr1, array3d * arr2, array3d *result);
void AApow(array3d * arr1, double M, array3d *result);
void Asum2(array3d * arr1, array3d *result);
void Asum3(array3d * arr1, matrix *result);
void Msum1(matrix * arr1, row *result);
void MMproddot(matrix * arr1, matrix * arr2, double M, matrix *result);
void MMsum(matrix * arr1, matrix * arr2, matrix *result);
void AAprodscalar(array3d * arr1, double M, array3d *result);
void AAexpminus(array3d * arr1, array3d *result);
double Vprodsqrt( row* arr1);
void VVsum(row * arr1, row * arr2, row *result);
void ASprodminus( array3d* arr1, double s, double m, array3d *result);
void VVproddot(row * arr1, row* arr2, double M, row *result);

void write_txtf_M(matrix *V, char opfile[]);
void write_txtf(vector<double> * V, char opfile[]);

void MtoA (matrix *M, array3d *A, int d2);
void read_txtf(matrix *data, char filename[]);
void tM (matrix *Min, matrix *Mout);


//#include "vargplvmheader.h"
/*void repmatRow(row * vec, int x, int y, int z, array3d *arr)
{
    matrix mat(x, (*vec));
    if (y == 1)
    {
        for (int i = 0; i < x; i++)
        {
            mat[i] = (*vec);
        }
    }
    (*arr).resize(z, mat);
    for (int i = 0; i < z; i++)
    {
        (*arr)[i] = mat;
    }
}
void repmatArr(array3d * arr, int x, int y, int z, array3d *tmp)
{
    int n, m;
    n = (*arr).size();
    m = (*arr)[0].size();

    row vec(y);
    matrix mat(m,vec);
    (*tmp).resize(n, mat);

    if (x == 1 && z == 1)
    {
        for (int ii = 0; ii < n; ii++)
        {
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < y; i++)
                {
                    (*tmp)[ii][j][i] = (*arr)[ii][j][0];
                }
            }

        }
    }
}
void  AAminus(array3d * arr1, array3d * arr2, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
            {
                (*result)[i][j][k] = (*arr1)[i][j][k] - (*arr2)[i][j][k];
            }
        }
    }
}
void  AAdiv(array3d * arr1, array3d * arr2, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = (*arr1)[i][j][k] / (*arr2)[i][j][k];
        }
    }
}
void  AAprod(array3d * arr1, array3d * arr2, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = (*arr1)[i][j][k] * (*arr2)[i][j][k];
        }
    }
}
void  AAsum(array3d * arr1, array3d * arr2, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = (*arr1)[i][j][k] + (*arr2)[i][j][k];
        }
    }
}
void  AApow(array3d * arr1, double M, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = pow((*arr1)[i][j][k], M);
        }
    }
}
void  Asum2(array3d * arr1, array3d *result)
{
    row vec(1, 0.0);
    matrix mat((*arr1)[0].size(), vec);
    (*result).resize((*arr1).size(), mat);
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
            {
                if (k == 0)
                    (*result)[i][j][0] = (*arr1)[i][j][k];
                else
                    (*result)[i][j][0] = (*result)[i][j][0] + (*arr1)[i][j][k];
            }

        }
    }
}
void Asum3(array3d * arr1, matrix *result)
{
    for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
    {
        for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
        {
            for(unsigned int i = 0; i < (*arr1).size(); i++)
            {
                if (i == 0)
                    (*result)[j][k] = (*arr1)[i][j][k];
                else
                    (*result)[j][k] = (*result)[j][k] + (*arr1)[i][j][k];
            }

        }
    }
}
void Msum1(matrix * arr1, row *result)
{
    for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
    {
        for(unsigned int k = 0; k <(*arr1).size(); k++)
        {
            if (k == 0)
                (*result)[j] = (*arr1)[k][j];
            else
                (*result)[j] = (*result)[j] + (*arr1)[k][j];
        }
    }
}
void MMproddot(matrix * arr1, matrix * arr2, double M, matrix *result)
{
    for(unsigned int j = 0; j <(*arr1).size(); j++)
    {
        for(unsigned int k = 0; k <(*arr1)[0].size(); k++)
        {
            (*result)[j][k] = (*arr1)[j][k] * (*arr2)[j][k] * M;
        }
    }
}
void MMsum(matrix * arr1, matrix * arr2, matrix *result)
{
    for(unsigned int j = 0; j <(*arr1).size(); j++)
    {
        for(unsigned int k = 0; k <(*arr1)[0].size(); k++)
        {
            (*result)[j][k] = (*arr1)[j][k] +(*arr2)[j][k];
        }
    }
}
void AAprodscalar(array3d * arr1, double M, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = (*arr1)[i][j][k]*M;
        }
    }
}
void AAexpminus(array3d * arr1, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = exp(-(*arr1)[i][j][k]);
        }
    }
}
double Vprodsqrt( row* arr1)
{
    double result = 1.0;
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        result = result * sqrt((*arr1)[i]);
    }
    return result;
}
void  VVsum(row * arr1, row * arr2, row *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        (*result)[i] = (*arr1)[i] + (*arr2)[i];
    }
}
void ASprodminus( array3d* arr1, double s, double m, array3d *result)
{
    for(unsigned int i = 0; i < (*arr1).size(); i++)
    {
        for(unsigned int j = 0; j <(*arr1)[0].size(); j++)
        {
            for(unsigned int k = 0; k <(*arr1)[0][0].size(); k++)
                (*result)[i][j][k] = 2*((*arr1)[i][j][k]) - 1;
        }
    }
}

void  VVproddot(row * arr1, row* arr2, double M, row *result)
{
    for(unsigned int j = 0; j <(*arr1).size(); j++)
    {
        (*result)[j] = (*arr1)[j] * (*arr2)[j] * M;
    }
}

//#include "vargplvmheader.h"

void read_txtf(matrix *data, char filename[])
{
	char buffer[1000];
    char *p; // p is the pointer point to the first character of buffer 

	int j=0;// i and j are row and column indeces of c, which are start from 0 to 2503 and 99, respectively
	int count=0; // count for the number of ' '
    int col = 0;
    FILE *fp=fopen(filename,"r");
    
    if( !fp)
    {
    cout<<"Can't open file "<< filename<< ", exiting ...\n";
    cin.get();
    exit(1);
    }
    
    
    fgets(buffer, 1000,fp);
    p = buffer;
    while (*p!='\n')
    {
     p++;
     if (*p == '\t')
     col++;
    }
    (*data).resize(col+1);
    while( 1 )
    {
    char buffer[1000] = {'\0'};
    char buffer1[100] = {'\0'};
    fgets(buffer, 1000, fp);
      p = buffer; // the pointer 'p' point to the address of the first character of buffer
      if(feof(fp))
      break;
      while (*p != '\n')
      {
            if(*p == '\t')// space or not?
    		{
            buffer1[j]='\0';
            (*data)[count].push_back(atof(buffer1)); // convert string to float
   			count++;
   			j = 0;
   			p++;
    		}
            else 
            {
            buffer1[j]= *p; // to be stored in column 1
	        j++;
      	    p++;   
            }
    }
    if(*p == '\n')
    	    {		
            buffer1[j] = '\0';
            (*data)[count].push_back(atof(buffer1)); 
            j = 0;
    	    } 
    count = 0;
    j=0;
    }
    fclose(fp);
}

void MtoA (matrix *M, array3d *A, int d2)
{
     int d3 = (*M)[0].size()/d2;
     row vec(d2, 0.0);
     matrix tmp((*M).size(), vec);
     for (int i = 0; i < d3; i++)
     {
         for (int j = 0; j < d2; j++)
         {
             for (unsigned int x = 0; x < (*M).size(); x++)
             tmp[x][j] = (*M)[x][j+i*d2];
         }
         (*A)[i] = tmp;
     }
 }
 
void write_txtf(vector<double> * V, char opfile[])
{
    FILE *out;

    out = fopen(opfile,"w");
    for (unsigned int t = 0; t<(*V).size(); t++)
    {
        fprintf(out, "%.17g", (*V)[t]);
        if (t<=(*V).size()-1)
            fprintf(out, "\n");
    }
    fclose(out);
} 
void write_txtf_M(matrix *V, char opfile[])
{
    FILE *out;

    out = fopen(opfile,"w");
    if ((*V).size()>0)
    {
        for (unsigned int t = 0; t<(*V)[0].size(); t++)
        {
            for(unsigned int i = 0; i<(*V).size(); i++)
            {
    	        fprintf(out, "%.17g", (*V)[i][t]);
                if (i < (*V).size()-1)
    	           fprintf(out, "\t");
            }
            if (t <=(*V)[0].size()-1)
                fprintf(out, "\n");
        }
    }
    fclose(out);
}
*/
