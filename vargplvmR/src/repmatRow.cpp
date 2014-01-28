#include "vargplvmheader.h"
void repmatRow(row * vec, int x, int y, int z, array3d *arr)
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

