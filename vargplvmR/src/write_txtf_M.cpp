#include "vargplvmheader.h"

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


