#include "vargplvmheader.h"

void MtoA (matrix *M, array3d *A, int d2)
{
     int d3 = (*M).size()/d2;
     row vec(d2, 0.0);
     matrix tmp((*M)[0].size(), vec);
     for (int i = 0; i < d3; i++)
     {
         for (int j = 0; j < d2; j++)
         {
             for (unsigned int x = 0; x < (*M)[0].size(); x++)
             tmp[x][j] = (*M)[j+i*d2][x];
         }
         (*A)[i] = tmp;
     }
 }
