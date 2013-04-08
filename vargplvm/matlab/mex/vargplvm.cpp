/* INSTRUCTIONS:
 * Compiling to handle exceptions (that's for typedef): for the watcom compiler
 * windows: mex COMPFLAGS="$COMPFLAGS -xs" vargplvm.cpp
 */


#include "vargplvmheader.h"
#include "mex.h"
/*
 * #include <R.h>
 * #include <Rdefines.h>
 * #include <Rinternals.h>
 * #include <Rmath.h>
 */

int DEBUGno =0; ////// DEBUG (init)


//row-wise
/*
mxArray*  serializeMatrix(matrix mat) {
    std::vector<double> v;
    for (int i=0; i < mat.size(); i++) {
        for (int j=0; j < mat[0].size(); j++) {
            v.push_back(mat[i][j]);
        }
    }
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    
    return mx;
}
*/


//column-wise
mxArray*  serializeMatrix(matrix mat) {
    std::vector<double> v;
    for (unsigned int j=0; j < mat[0].size(); j++) {
        for (unsigned int i=0; i < mat.size(); i++) {
            v.push_back(mat[i][j]);
        }
    }
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
	// memory copy so that matlab doesnt free the memory after this point
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    
    return mx;
}


mxArray*  copyRow(row r) {
    std::vector<double> v;
    for (unsigned int i=0; i < r.size(); i++) {
            v.push_back(r[i]);
        }
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    
    return mx;
}


// Keep in mind that MATLAB is column major, C++ is row major
matrix reshapeMatrix(double* A, int N, int Q) {
    int glob =0;
    row vecZero(Q,0.0);
    matrix mat(N,vecZero);
    
    for (unsigned int i=0; i<N; i++){
        for (unsigned int j=0; j<Q; j++) {
            mat[i][j] = A[glob++];
        }
    }
    return mat;
}

array3d reshapeMatrix(double* A, int N, int Q, int M) {
    int glob =0;
        row vecZero(Q,0.0);
    matrix matZero(N,vecZero);
    array3d arr3d(M,matZero);
    
    for (unsigned int i=0; i<M ; i++) {
        for (unsigned int j=0; j<N; j++){
            for (unsigned int k=0; k<Q; k++) {
                arr3d[i][j][k] = A[glob++];
            }
        }
    }
    return arr3d;
}

void printMatrix(matrix mat) {
    for (unsigned int i=0; i < mat.size(); i++) {
        for (unsigned int j=0; j < mat[0].size(); j++) {
            printf("%f ",mat[i][j]);
        }
        printf("\n");
    }
}

void printMatrix(array3d ar3d) {
    
    for (unsigned int i=0; i<ar3d.size(); i++){
        for (unsigned int j=0; j<ar3d[0].size(); j++) {
            for (unsigned int k=0; k<ar3d[0][0].size(); k++) {
                printf("%f ",ar3d[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


void DEBUGwrite()
{
    FILE *out;
    
    out = fopen("DEBUG.txt","w");
    
    fprintf(out, "DEBUG%d\n",DEBUGno );
    DEBUGno++;
    fclose(out);
}

/**********************/

void read_txtf(matrix *data, char filename[])
{
    char buffer[1000000];
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
    
    fgets(buffer, 1000000,fp);
    p = buffer;
    while (*p!='\n')
    {
        p++;
        if (*p == '\t')
            col++;
    }
    (*data).resize(col+1);
    //cout<<"data size "<<(*data).size() << " ";
    while( 1 )
    {
        char buffer[1000000] = {'\0'};
        char buffer1[1000] = {'\0'};
        fgets(buffer, 1000000, fp);
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

//Matrix to array
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


// Transpose matrix
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




/****************************************/


//void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray *prhs[]){
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *MNQ, *A_vec, *covGrad_dim2, *means_arr, *covars_arr, *asPlus1_arr, *aDasPlus1_arr, *ZmZm_arr, *covGrad_arr;
    int M,N,Q, le, Zd2, cd2;
     
    
    /* Check for proper number of arguments */
    if (nrhs != 9) {
        mexErrMsgIdAndTxt( "MATLAB:vargplvm:invalidNumInputs",
                "Nine input arguments required (double* [M, N, Q], double* A_vec, double* covGrad_dim2, double* meansTransp_arr, double* covarsTransp_arr, double* asPlus1Transp_arr, double* aDasPlus1Transp_arr, double* ZmZm_arr, double* covGrad_arr).");
    } else if (nlhs > 4) {
        mexErrMsgIdAndTxt( "MATLAB:yprime:maxlhs",
                "Too many output arguments.");
    }
    
    
//mwSize i, n;
    MNQ = mxGetPr(prhs[0]);
    A_vec = mxGetPr(prhs[1]);
    covGrad_dim2 = mxGetPr(prhs[2]);
    means_arr = mxGetPr(prhs[3]);
    covars_arr = mxGetPr(prhs[4]);
    asPlus1_arr = mxGetPr(prhs[5]);
    aDasPlus1_arr = mxGetPr(prhs[6]);
    ZmZm_arr = mxGetPr(prhs[7]);
    covGrad_arr = mxGetPr(prhs[8]);
    
    M = (int) MNQ[0];
    N = (int) MNQ[1];
    Q = (int) MNQ[2];
    
//ZmZm_arr = mxGetPr(prhs[2]);
    Zd2 = Q; //(int) ZmZm_arr[0];
    cd2 = (int) covGrad_dim2[0];
    
    //mexPrintf("%f\n",A_vec[0]);///
    //mexPrintf("M=%d, N=%d, Q=%d\n",M,N,Q);///////
    //mexPrintf("Zd2=%d, cd2=%d\n",Zd2, cd2); ///////
    //mexPrintf("%d\n",nrhs);//////
    
    double *pA;
    pA = A_vec; // !!!!! pA = real(A_vec);
    //le = Q;
    
    //int M, N, Q, le, Zd2, cd2;
    ////bool isarray;
    //M = INTEGER_VALUE(M_int);
    //N = INTEGER_VALUE(N_int);
    //Q = INTEGER_VALUE(Q_int);
    //Zd2 = INTEGER_VALUE(ZmZm_arr);
    //cd2 = INTEGER_VALUE(covGrad_arr);
    
    //double *pA;
    //pA = REAL(A_vec);
    //le = length(A_vec);
    
    row A(Q, 0.0);
    for (int i = 0; i < Q; i++)
        A[i] = pA[i];
    row ve(Q);
    
    matrix means;
    means = reshapeMatrix(means_arr, N,Q);

    matrix covars;
    covars = reshapeMatrix(covars_arr, N,Q);
       
    matrix asPlus1;
    asPlus1=reshapeMatrix(asPlus1_arr, N,Q);

    matrix aDasPlus1;
    aDasPlus1=reshapeMatrix(aDasPlus1_arr, N,Q);

    //row vecZeroM(M*Q,0.0);
    matrix dataA;//(M,vecZeroM);
   // char Zm1[] = "ZmZm.txt";
    //read_txtf(&dataA, Zm1);
    dataA = reshapeMatrix(ZmZm_arr, M*Q,M);
    //printMatrix(dataA);////
    
    matrix data2(M,ve);
    array3d ZmZm(M,data2);
    MtoA (&dataA, &ZmZm, Zd2);
    // printf("ZmZm in cpp %u, %u, %u, %f\n", ZmZm.size(), ZmZm[0].size(), ZmZm[0][0].size(), ZmZm[0][0][0]);////////////
    
    int n1 = 0, m1 = 0,x1 = 0;
   
    matrix dataB;
    //char cov1[] = "covGrad.txt";
    //read_txtf(&dataB, cov1);
    dataB = reshapeMatrix(covGrad_arr, M*cd2, M);
    
    m1 = cd2;
    n1 = dataB.size();
    x1 = dataB[0].size()/m1;
    row ve3(m1,0);
    matrix data3(n1,ve3);
    array3d covGrad(x1,data3);
    MtoA (&dataB, &covGrad, cd2);
    
    // printf("covGrad in cpp %u, %u, %u, %.16f\n", covGrad.size(), covGrad[0].size(), covGrad[0][0].size(), covGrad[0][0][0]);//////////
    
    //getArray(covGrad_arr, &covGrad);
    
    row mu_n(Q);
    row s2_n(Q);
    row AS_n(Q);

    array3d MunZmZm(ZmZm);
    array3d MunZmZmA(ZmZm);
    array3d k2Kern_n(covGrad);
    array3d tmpA(ZmZm);
    array3d k2ncovG(ZmZm);
    
    array3d *pMunZmZm = &MunZmZm;
    array3d *pMunZmZmA = &MunZmZmA;
    array3d *pk2Kern_n = &k2Kern_n;
    array3d *ptmpA = &tmpA;
    array3d *pk2ncovG = &k2ncovG;

    array3d ar2k(covGrad);
    array3d *par2k = &ar2k;
    
    matrix tmpM(M,A);
    matrix *tmp = &tmpM;
    matrix gVarcovars(N,mu_n);
    matrix gVarmeans(gVarcovars);
    matrix partInd2(M,mu_n);
    
    row partA2(Q);
    matrix Amq(M, A);
    array3d ar(ZmZm);
    array3d *par = &ar;
// printf("par in cpp %u, %u, %u, %f\n", (*par).size(), (*par)[0].size(), (*par)[0][0].size(), (*par)[0][0][0]);

    array3d ar2(ZmZm);
    array3d *par2 = &ar2;
    
    array3d ar3(ZmZm);
    array3d *par3 = &ar3;
    
    array3d ar4(ZmZm);
    array3d *par4 = &ar4;
    array3d ar5(ZmZm);
    array3d *par5 = &ar5;
    
    row arV(Q);
    row *parV =&arV;
    row arV2(Q);
    row *parV2 = &arV2;
    
    matrix arM(M,A);
    matrix *parM = &arM;
    
    matrix arM2(M,A);
    matrix *parM2 = &arM2;
    

//printf("here \n");
    
    for (int n = 0; n < N; ++n)
    {
        
        mu_n = means[n];
        s2_n = covars[n];
        AS_n = asPlus1[n];
        //printf("here 20\n");
        repmatRow(&mu_n, M, 1, M, par);
        //printf("par in cpp %u, %u, %u, %.16f\n", (*par).size(), (*par)[0].size(), (*par)[0][0].size(), (*par)[0][0][0]);
        
        AAminus(par, &ZmZm, pMunZmZm);
        //printf("here 21\n");
        repmatRow(&AS_n, M, 1, M, par);
        AAdiv(pMunZmZm, par, pMunZmZmA);
        //printf("here 22\n");
        
        AApow(pMunZmZm, double(2), par2);
        repmatRow(&(aDasPlus1[n]), M, 1, M, par);
        AAprod(par2, par , par3);
        Asum2(par3,pk2Kern_n);
        //printf("here 23\n");
        
        AAexpminus(pk2Kern_n, par2k);
        AAprodscalar(par2k, 1/Vprodsqrt(&AS_n), pk2Kern_n);
        
        // derivatives wrt to variational means
        AAprod(pk2Kern_n,&covGrad, par2k);
        repmatArr(par2k, 1, Q, 1, pk2ncovG);
        
        AAprod(pMunZmZmA,pk2ncovG,ptmpA);
        Asum3(ptmpA, tmp);
        //printf("here 24\n");
        Msum1(tmp, parV2);
        VVproddot(&A, parV2, double(-2), parV);
        //printf("here 25\n");
        gVarmeans[n] = *parV;
        // derivatives wrt inducing inputs
        MMproddot(&Amq, tmp, double(1), parM);
        
        MMsum(&partInd2, parM, parM2);
        
        partInd2 = *parM2;
        // Derivative wrt input scales
        AAprod(pMunZmZmA, pMunZmZm, pMunZmZmA);
        //printf("here 26\n");
        repmatRow(&s2_n, M, 1, M, par);
        AAsum(pMunZmZmA, par, par2);
        repmatRow(&AS_n, M, 1, M, par3);
        AAprod(par2, pk2ncovG, par4);
        AAdiv(par4,par3 , par5);
        Asum3(par5, parM);
        Msum1(parM,parV);
        VVsum(&partA2, parV, parV2);
        partA2 = *parV2;
        //           printf("here 27\n");
        // derivatives wrt variational diagonal covariances
        repmatRow(&A, M, 1, M, par);
        AAprod(pMunZmZmA, par,pMunZmZmA);
        repmatRow(&(aDasPlus1[n]), M, 1, M, par);
        ASprodminus(pMunZmZmA, double (2), double(1), par4);
        AAprod(par4, pk2ncovG,par3);
        AAprod(par,    par3   ,par2);
        Asum3(par2,parM);
        Msum1(parM,parV);
        gVarcovars[n] = *parV;
        //printf("here 2\n");
    }

//printf("here 3\n");
    
    // write to text
    /*
    char pI2[] = "partInd2.txt";
    write_txtf_M(&partInd2, pI2);
    char pA2[] = "partA2.txt";
    write_txtf(&partA2, pA2);
    char gVm[] = "gVarmeans.txt";
    write_txtf_M(&gVarmeans, gVm);
    char gVc[] = "gVarcovars.txt";
    write_txtf_M(&gVarcovars, gVc);
    */
    
    ///////////
    plhs[0] = serializeMatrix(partInd2);
    plhs[1] = copyRow(partA2);
    plhs[2] = serializeMatrix(gVarmeans);
    plhs[3] = serializeMatrix(gVarcovars);
    //return 1;
}


