#include "vargplvmheader.h"

void tM (matrix *Min, matrix *Mout)
{
    row A((*Min).size(), 0.0);
    (*Mout).resize((*Min)[0].size(), A);
     for (int i = 0; i < (*Min).size() ; i++)
     {
         for (int j = 0; j < (*Min)[0].size(); j++)
         {
             (*Mout)[j][i] = (*Min)[i][j];
         }
     }
 }
