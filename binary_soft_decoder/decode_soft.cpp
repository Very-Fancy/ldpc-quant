#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>
#include <list>
#include <iostream>
#include <mex.h>
#include <matrix.h>
#include <conio.h> 
#include <fstream>
#include "Matrix.h"
#include "common.h"
#include "LDPC.h"
#include "ObjectHandle.h"

#define ABS(A) (((A) >= 0) ? (A) : -(A))
#define HARD(A) (((A) < 0)?1:0)
#define OFFSET(A, B) (((A) > (B)) ? ((A)-(B)) : 0)
#define INFTY 1000000

using namespace std;

int clip(double val, int qn)
{
	int res = (int)round(val);
	if (res < -(int)pow(2, qn - 1) )
	{
		return -(int)(pow(2,(double)qn - 1));
	}
	else
	{
		if (res > (int)pow(2,(double)qn - 1) - 1)
		{
			return  (int)pow(2,(double)qn - 1) - 1;
		}
		else
		{
			return (int)res;
		}
	}
}

bool MinSum(LDPC& ldpc, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<double>& scale_array, vector<double>& offset_array)
{
    // Auxiliary matrices
    Matrix<double>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<double>  Q_msgs(ldpc.m, ldpc.rmax); // messages from variable to check nodes

    // Initialization
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);

        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = in_llr[i];
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            double first_min = INFTY;
            double second_min = INFTY;
            int min_index = 0;
            int sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                double temp = fabs(Q_msgs(j, k));
                if (temp < first_min)
                {
                    second_min = first_min;
                    first_min = temp;
                    min_index = k;
                }
                else if (temp < second_min)
                {
                    second_min = temp;
                }

                if (Q_msgs(j, k) < 0)
                {
                    sum_sign ^= 1;
                }
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                int sign = sum_sign;
                if (Q_msgs(j, k) < 0)
                {
                    sign ^= 1;
                }
                if (k == min_index)
                {
                    R_msgs(j,k) = second_min;
                }
                else
                {
                    R_msgs(j,k) = first_min;
                }
                double scale = scale_array[ldpc.row_col(j, k)];
                double offset = offset_array[ldpc.row_col(j, k)];
                R_msgs(j,k) = (1-2*sign)*scale*OFFSET(R_msgs(j,k), offset);
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] += R_msgs(ldpc.msgs_col(i,k));
            }

            y[i] = HARD(out_llr[i]);

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                Q_msgs(ldpc.msgs_col(i,k)) = out_llr[i] - R_msgs(ldpc.msgs_col(i,k));
            }
        }
    }

    if (number_of_iter)
    {
        *number_of_iter = max_iter;
    }
    return 1;
}

bool MinSumQ(LDPC& ldpc, double qn, vector<double>& in_llr, int max_iter, vector<double>& out_llr, vector<int>& y, int* number_of_iter, vector<double>& scale_array, vector<double>& offset_array)
{
   
    Matrix<int>  R_msgs(ldpc.m, ldpc.rmax); // messages from check to variable nodes
    Matrix<int>  Q_msgs(ldpc.m, ldpc.rmax); // messages from variable to check nodes
  
    for (int i = 0; i < ldpc.n; ++i)
    {
        out_llr[i] = in_llr[i];
        y[i] = HARD(out_llr[i]);

        for (int j = 0; j < ldpc.col_weight[i]; ++j)
        {
            Q_msgs(ldpc.msgs_col(i,j)) = clip(in_llr[i], qn);
        }
    }

    for (int loop = 0; loop < max_iter; ++loop)
    {
        // Update R
        for (int j = 0; j < ldpc.m; j++)
        {
            double first_min = INFTY;
            double second_min = INFTY;
            int min_index = 0;
            int sum_sign = 0;

            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                double temp = fabs(Q_msgs(j, k));
                if (temp < first_min)
                {
                    second_min = first_min;
                    first_min = temp;
                    min_index = k;
                }
                else if (temp < second_min)
                {
                    second_min = temp;
                }

                if (Q_msgs(j, k) < 0)
                {
                    sum_sign ^= 1;
                }
            }
            for (int k = 0; k < ldpc.row_weight[j]; k++)
            {
                int sign = sum_sign;
                if (Q_msgs(j, k) < 0)
                {
                    sign ^= 1;
                }
                if (k == min_index)
                {
                    R_msgs(j,k) = (int)second_min;
                }
                else
                {
                    R_msgs(j,k) = (int)first_min;
                }
             
                double scale = scale_array[ldpc.row_col(j, k)];
                double offset = offset_array[ldpc.row_col(j, k)];
                R_msgs(j,k) = (int)clip(round((1-2*sign)*scale*OFFSET(R_msgs(j,k), offset)), qn);
            }
        }
        // Update Q
        for (int i = 0; i < ldpc.n; i++)
        {
            out_llr[i] = in_llr[i];

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                out_llr[i] = out_llr[i] + R_msgs(ldpc.msgs_col(i,k));
            }

            y[i] = HARD(out_llr[i]);

            for (int k = 0; k < ldpc.col_weight[i]; k++)
            {
                double tmp = in_llr[i];
                for (int ii = 0; ii < ldpc.col_weight[i]; ii++)
                {
                    if (ii != k)
                    {
                        tmp += R_msgs(ldpc.msgs_col(i,ii) );
                    }
                 }
                Q_msgs(ldpc.msgs_col(i,k)) = clip(tmp,qn);
            }  
        }
    }

    if (number_of_iter)
    {
        
        *number_of_iter = max_iter;
    }
    return 1;
}

//=======================================================================================================================================
//================================================================== Mex function =======================================================
//=======================================================================================================================================


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int command = (int)*(mxGetPr(prhs[0]));
    double* input = 0;

    switch(command)
    {
        case -2:
        {
            char buf[256];
            mxGetString(prhs[1], buf, 256); 
             
            LDPC* pLdpc = new LDPC();
            pLdpc->init(buf);
            
            if (nlhs > 0)
            {
	            plhs[0] = create_handle(pLdpc);
            }
            if (nlhs > 1)
            {
                plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[1]))) = pLdpc->n;
            }
            if (nlhs > 2)
            {	
                plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
                *((double*)(mxGetPr(plhs[2]))) = pLdpc->m;
            }
            break;   
        }
        case -1:
        {
            // LDPC object
            destroy_object<LDPC>(prhs[1]);
            break;   
        }
        /* Sum-Product */
       
        /* Min-Sum */
        case 0:
        {
          // LDPC object
            LDPC& ldpc = get_object<LDPC>(prhs[1]);

            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<double> scale_array(ldpc.n);
            vector<double> offset_array(ldpc.n);
            vector<int> y(ldpc.n);

            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // scale array
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                scale_array[i] = input[i];
            }

            // offset array
            input = mxGetPr(prhs[5]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                offset_array[i] = input[i];
            }

            bool result = true;
            result = MinSum(ldpc, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, scale_array, offset_array);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = out_llr[i];
            }
            break;
        }
        default:
        {
            int qn = command;
              LDPC& ldpc = get_object<LDPC>(prhs[1]);
     //       double qn = 6;
            vector<double> in_llr(ldpc.n);
            vector<double> out_llr(ldpc.n);
            vector<double> scale_array(ldpc.n);
            vector<double> offset_array(ldpc.n);
            vector<int> y(ldpc.n);
           // int qn = (int)*(mxGetPr(prhs[2]));
           
            // in LLR
            input = mxGetPr(prhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                in_llr[i] = input[i];
            }

            // iterations
            int iMaxNumberOfIterations = (int)*(mxGetPr(prhs[3]));
            int number_of_iter = 0;

            // scale array
            input = mxGetPr(prhs[4]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                scale_array[i] = input[i];
            }

            // offset array
            input = mxGetPr(prhs[5]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                offset_array[i] = input[i];
            }

            bool result = true;
            result = MinSumQ(ldpc, qn, in_llr, iMaxNumberOfIterations, out_llr, y, &number_of_iter, scale_array, offset_array);

            /* denial flag */
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[0]))) = result;

            /* number of iterations */
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            *((double*)(mxGetPr(plhs[1]))) = number_of_iter;

            double* data = NULL;

            /* hard desicion */
            plhs[2] = mxCreateDoubleMatrix(1, ldpc.n, mxREAL);
            data = mxGetPr(plhs[2]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = y[i];
            }

            /* out LLR */
            plhs[3] = mxCreateDoubleMatrix (1,ldpc.n ,mxREAL );
            data = mxGetPr(plhs[3]);
            for(int i = 0; i < ldpc.n; ++i)
            {
                data[i] = clip(out_llr[i],qn);
            }
            break;
            break;
        }
               
       
            // LDPC object
           
        
       
    }

}