#include "vargplvmheader.h"
extern "C" {

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

    //void getArray(SEXP Rvec, array3d * arr);
//    void getMatrix(SEXP Rvec, matrix * data);
   // SEXP covertMatrix(matrix *inmat);
    //SEXP covertVector(row *inmat);
    SEXP vargplvm(SEXP M_int, SEXP N_int, SEXP Q_int, SEXP A_vec,
                  //SEXP means_mat, SEXP covars_mat, SEXP asPlus1_mat, SEXP aDasPlus1_mat,
                  SEXP ZmZm_arr, SEXP covGrad_arr);

    SEXP vargplvm(SEXP M_int, SEXP N_int, SEXP Q_int, SEXP A_vec,
                 // SEXP means_mat, SEXP covars_mat, SEXP asPlus1_mat, SEXP aDasPlus1_mat,
                  SEXP ZmZm_arr, SEXP covGrad_arr)
    {
        int M, N, Q, le, Zd2, cd2;
        //bool isarray;
        M = INTEGER_VALUE(M_int);
        N = INTEGER_VALUE(N_int);
        Q = INTEGER_VALUE(Q_int);
        Zd2 = INTEGER_VALUE(ZmZm_arr);
        cd2 = INTEGER_VALUE(covGrad_arr);

        double *pA;
        pA = REAL(A_vec);
        le = length(A_vec);
        row A(le, 0.0);
        for (int i = 0; i < le; i++)
            A[i] = pA[i];
        row ve(Q);


        matrix tmeans, means;
        char mea[] = "means.txt";
        read_txtf(&tmeans, mea);
        tM(&tmeans, &means);
       // printf("means in cpp %lu, %u, %f\n", means.size(), means[0].size(), means[0][0]);
        //getMatrix(means_mat, &means);
        matrix covars, tcovars;
        char cov[] = "covars.txt";
        read_txtf(&tcovars, cov);
        tM(&tcovars, &covars);
         //       printf("covars in cpp %lu, %u, %f\n", covars.size(), covars[0].size(), covars[0][0]);
        //getMatrix(covars_mat, &covars);
        matrix asPlus1, tasPlus1;
        char asP[] = "asPlus1.txt";
        read_txtf(&tasPlus1, asP);
        tM(&tasPlus1, &asPlus1);
         //       printf("asPlus1 in cpp %lu, %u, %f\n", asPlus1.size(), asPlus1[0].size(), asPlus1[0][0]);
        //getMatrix(asPlus1_mat, &asPlus1);
        matrix aDasPlus1, taDasPlus1;
        char aDa[] = "aDasPlus1.txt";
        read_txtf(&taDasPlus1, aDa);
        tM(&taDasPlus1, &aDasPlus1);
         //       printf("aDasPlus1 in cpp %lu, %u, %f\n", aDasPlus1.size(), aDasPlus1[0].size(), aDasPlus1[0][0]);
        //getMatrix(aDasPlus1_mat, &aDasPlus1);

        matrix dataA;
        char Zm1[] = "ZmZm.txt";
        read_txtf(&dataA, Zm1);
        //printf("dataA in cpp %lu, %u, %f\n", dataA.size(), dataA[0].size(), dataA[0][0]);

        matrix data2(M,ve);
        array3d ZmZm(M,data2);
        MtoA (&dataA, &ZmZm, Zd2);
        //printf("ZmZm in cpp %u, %u, %u, %f\n", ZmZm.size(), ZmZm[0].size(), ZmZm[0][0].size(), ZmZm[0][0][0]);
        //getArray(ZmZm_arr, &ZmZm);

        //SEXP Rdim = getAttrib(covGrad_arr, R_DimSymbol);
        int n1 = 0, m1 = 0,x1 = 0;
        //n1 = INTEGER(Rdim)[0];
        //m1 = INTEGER(Rdim)[1];
        //x1 = INTEGER(Rdim)[2];

        matrix dataB;
        char cov1[] = "covGrad.txt";
        read_txtf(&dataB, cov1);
        //printf("dataB in cpp %lu, %u, %.16f\n", dataB.size(), dataB[0].size(), dataB[0][0]);
        m1 = cd2;
        n1 = dataB.size();
        x1 = dataB[0].size()/m1;
        row ve3(m1,0);
        matrix data3(n1,ve3);
        array3d covGrad(x1,data3);
        MtoA (&dataB, &covGrad, cd2);
        //printf("covGrad in cpp %u, %u, %u, %.16f\n", covGrad.size(), covGrad[0].size(), covGrad[0][0].size(), covGrad[0][0][0]);

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
        char pI2[] = "partInd2.txt";
        write_txtf_M(&partInd2, pI2);
        char pA2[] = "partA2.txt";
        write_txtf(&partA2, pA2);
        char gVm[] = "gVarmeans.txt";
        write_txtf_M(&gVarmeans, gVm);
        char gVc[] = "gVarcovars.txt";
        write_txtf_M(&gVarcovars, gVc);
//printf("here 4\n");
        // return all results in SEXP list
        /*int nProtected = 0;
        SEXP partInd2_mat, partA2_mat, gVarmeans_mat, gVarcovars_mat;
        SEXP list, names;

        partInd2_mat = covertMatrix(&partInd2);

        partA2_mat = covertVector(&partA2);

        gVarmeans_mat = covertMatrix(&gVarmeans);

        gVarcovars_mat = covertMatrix(&gVarcovars);

        PROTECT(list = allocVector(VECSXP, 4));
        ++nProtected;
        SET_VECTOR_ELT(list, 0, partInd2_mat);
        SET_VECTOR_ELT(list, 1, partA2_mat);
        SET_VECTOR_ELT(list, 2, gVarmeans_mat);
        SET_VECTOR_ELT(list, 3, gVarcovars_mat);
        PROTECT(names = allocVector(STRSXP, 4));
        ++nProtected;
        SET_STRING_ELT(names, 0, mkChar("partInd2"));
        SET_STRING_ELT(names, 1, mkChar("partA2"));
        SET_STRING_ELT(names, 2, mkChar("gVarmeans"));
        SET_STRING_ELT(names, 3, mkChar("gVarcovars"));
        setAttrib(list, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return list;*/
        return R_NilValue;
    }
   /* SEXP covertMatrix(matrix *inmat)
    {
        int nProtected = 0, nrows, ncols, i, j;
        SEXP dim, outmat;
        nrows = (*inmat).size();
        ncols = (*inmat)[0].size();
        PROTECT(dim = allocVector(INTSXP, 2));
        ++nProtected;
        PROTECT(outmat = allocVector(REALSXP, nrows * ncols));
        ++nProtected;

        INTEGER(dim)[0] = nrows;
        INTEGER(dim)[1] = ncols;

        for(i = 0; i < nrows; ++i)
        {
            for(j = 0; j < ncols; ++j)
            {
                REAL(outmat)[i + j*nrows] = (*inmat)[i][j];
            }
        }
        setAttrib(outmat, R_DimSymbol, dim);
        UNPROTECT(nProtected);
        return outmat;
    }

    SEXP covertVector(row *inmat)
    {
        int nProtected = 0, nrows, i;
        SEXP outmat;
        nrows = (*inmat).size();

        PROTECT(outmat = allocVector(REALSXP, nrows));
        ++nProtected;

        for(i = 0; i < nrows; ++i)
        {
            REAL(outmat)[i] = (*inmat)[i];
        }
        UNPROTECT(nProtected);
        return outmat;
    }
    void getArray(SEXP Rvec, array3d * arr)
    {
        int i, j, y = 0, n = 0, m = 0,x = 1, c = 0;
        double *vec;
        vec = REAL(Rvec);
        SEXP Rdim = getAttrib(Rvec, R_DimSymbol);
        int le = length(Rdim);

        n = INTEGER(Rdim)[0];
        m = INTEGER(Rdim)[1];

        x = INTEGER(Rdim)[2];

        for(y = 0; y <x; y++)
        {
            for(j = 0; j <m; j++)
            {
                for (i = 0; i < n; i++)
                {
                    (*arr)[y][i][j] = vec[c];
                    c++;
                }
            }
        }
    }
    void getMatrix(SEXP Rvec, matrix * data)
    {
        int i, j, n = 0, m = 0, c = 0;
        double *vec;
        vec = REAL(Rvec);
        SEXP Rdim = getAttrib(Rvec, R_DimSymbol);
        n = INTEGER(Rdim)[0];
        m = INTEGER(Rdim)[1];

        for(j = 0; j <m; j++)
        {
            for (i = 0; i < n; i++)
            {
                (*data)[i][j] = vec[c];
                c++;
            }
        }
    }*/
}
