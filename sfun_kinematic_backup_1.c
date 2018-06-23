/*
 * sfuntmpl_basic.c: Basic 'C' template for a level 2 S-function.
 *
 *  -------------------------------------------------------------------------
 *  | See matlabroot/simulink/src/sfuntmpl_doc.c for a more detailed template |
 *  -------------------------------------------------------------------------
 *
 * Copyright 1990-2002 The MathWorks, Inc.
 * $Revision: 1.27.4.2 $
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfuntmpl_basic with the name of your S-function).
 *A[0][0] = 2;A[0][1] = 3;A[0][2] = 2;
    A[1][0] = 3;A[1][1] = 1;A[1][2] = 3;
    B[0][0] = 2;B[0][1] = 3;
    B[1][0] = 3;B[1][1] = 1;
    B[2][0] = 0;B[2][1] = 1;
    cp_D[0][0] = 2;cp_D[0][1] = 3;
    cp_D[1][0] = 3;cp_D[1][1] = 1;
    MATRIX_Mul(A,B, D, 2, 3, 2);
    printf("D:%4f,%4f\n,%4f,%4f\n",D[0][0],D[0][1],D[1][0],D[1][1]);
    MATRIX_Add(D,cp_D,va_D,2,2);
    printf("Add:%4f,%4f\n,%4f,%4f\n",va_D[0][0],va_D[0][1],va_D[1][0],va_D[1][1]);
    MATRIX_Sub(D,cp_D,va_D,2,2);
    printf("Sub:%4f,%4f\n,%4f,%4f\n",va_D[0][0],va_D[0][1],va_D[1][0],va_D[1][1]);
    MATRIX_Tran(A,B,2,3);
     printf("B:%4f,%4f\n,%4f,%4f\n,%4f,%4f\n",B[0][0],B[0][1],B[1][0],B[1][1],B[2][0],B[2][1]);
    MATRIX_Pinv(A,B,2,3);
    printf("Pinv:%4f,%4f\n,%4f,%4f\n,%4f,%4f\n",B[0][0],B[0][1],B[1][0],B[1][1],B[2][0],B[2][1]);
 */

#define S_FUNCTION_NAME  sfun_kinematic
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include "mex.h" 
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

// 定义0矩阵， n*m
void MATRIX_SetZero(double *A, int n,int m)
{
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < m; j++){
            A[i*m+j] = 0.0;
        }
}
// 定义单位阵， n*n
void MATRIX_SetUnit(double *A, int n)
{
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++){
            if( i == j)
            {
                A[i*n+j] = 1.0;
            }
            else
            {
                A[i*n+j] = 0.0;
            }
        }
}
// 定义矩阵相乘 C = A * B,C为n*k，A为n*m，B为m*k
void MATRIX_Mul(double *A, double *B, double *C,int n,int m,int k)
{
    int i,j,p;
    for(i = 0; i < n; i++)
        for(j = 0; j < k; j++){
                C[i*k+j] = 0.0;
            for(p = 0; p < m; p++){
                C[i*k+j] = C[i*k+j] + (A[i*m+p]) * (B[p*k+j]);
        
            }
        }
}
void MATRIX_Add(double *A, double *B, double *C,int n,int m)//C = A + B
{
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < m; j++){
            C[i*m+j] = A[i*m+j] + B[i*m+j];
        }
}
void MATRIX_Sub(double *A, double *B, double *C,int n,int m)//C = A - B
{
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < m; j++){
            C[i*m+j] = A[i*m+j] - B[i*m+j];
        }
}
// 定义矩阵转置
void MATRIX_Tran(double *A, double *C,int n,int m)//C=A^T
{
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < m; j++){
            C[i+j*n] = A[i*m+j];
        }
}

// 方阵求逆
int MATRIX_Inv(double *C, double *IC, int n)//IC = inv(C)
{
	int i, j, k, l;
	
	/* 单位阵*/
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++) 
		{	
			*(IC+i*n+j) = 0.0;
		}
		*(IC+i*n+i) = 1.0;
	}
	
	/* 化上三角阵*/
	for (j=0; j<n; j++)
	{	
		if(fabs(*(C+j*n+j))>1e-15) /* C[j][j]不等于0*/
		{
			/* IC阵的第j行除以C[j][j]*/
			for(k=0; k<n; k++)
			{
				*(IC+j*n+k) /= *(C+j*n+j);
			}
			/* C阵的第j行除以C[j][j]*/
			for(k=n-1; k>=j; k--)
			{
				*(C+j*n+k) /= *(C+j*n+j);
			}
			
			for(i=j+1; i<n; i++)
			{
				/* IC阵的第i行 - C[i][j]*IC阵的第j行*/
				for(k=0; k<n; k++)
				{
					*(IC+i*n+k) -= *(C+i*n+j) * *(IC+j*n+k);
				}
				/* C阵的第i行 - C[i][j]*C阵的第j行*/
				for (k=n-1; k>=j; k--)
				{
					*(C+i*n+k) -= *(C+i*n+j) * *(C+j*n+k);
				}
			}
		}
		else if (j<n-1)
		{
			
			for(l=j+1; l<n; l++)
			{
				/* 若C阵第j行后的C[l][j]不等于0，第j行加上第l行*/
				if (fabs(*(C+l*n+j)) > 1e-15)
				{
					for (k=0; k<n; k++)
					{
						*(IC+j*n+k) += *(IC+l*n+k);
					}
					for (k=n-1; k>=j; k--)
					{
						*(C+j*n+k) += *(C+l*n+k);
					}
					/* IC阵的第j行除以C[j][j]*/
					for (k=0; k<n; k++)
					{
						*(IC+j*n+k) /= *(C+j*n+j);
					}
					/* C阵的第j行除以C[j][j]*/
					for (k=n-1; k>=j; k--)
					{
						*(C+j*n+k) /= *(C+j*n+j);
					}
					/* 第i行 - C[i][j]*第j行*/
					for (i=j+1; i<n; i++)
					{
						for (k=0; k<n; k++)
						{
							*(IC+i*n+k) -= *(C+i*n+j) * *(IC+j*n+k);
						}
						for (k=n-1; k>=j; k--)
						{
							*(C+i*n+k) -= *(C+i*n+j) * *(C+j*n+k);
						}
					}
					break;
				}
			}
			
			if (l == n)  /* C[l][j] 全等于0*/
			{
				return (-1);   /*  矩阵的行列式为零，不可求逆*/
			}
		}
		else  /* C[n][n]等于0*/
		{
			return (-1);    /*  矩阵的行列式为零，不可求逆*/
		}
	}
	/* 化成单位阵*/
	for(j=n-1; j>=1; j--)
	{
		for(i=j-1; i>=0; i--)
		{
			for(k=0; k<n; k++)
			{
				*(IC+i*n+k) -= *(C+i*n+j) * *(IC+j*n+k);
			}
			*(C+i*n+j) = 0;
		}
	}
	
	return (1);
}
// 求广义逆 C = A^T*inv(A*A^T)
void MATRIX_Pinv(double *A, double *C, int m, int n){
	double transA[10][10];
	double TempA[10][10];
	double invA[10][10];
	double TempB[10][10];
	double invB[10][10];
	MATRIX_Tran(A, transA, m, n);
   // printf("transA:%4f,%4f\n,%4f,%4f\n,%4f,%4f\n",transA[0][0],transA[0][1],transA[1][0],transA[1][1],transA[2][0],transA[2][1]);
	if(m < n){
		MATRIX_Mul(A, transA, TempB, m, n, m);
       // printf("TempB:%4f,%4f,%4f\n,%4f,%4f,%4f\n,%4f,%4f,%4f\n",TempB[0][0],TempB[0][1],TempB[0][2],TempB[1][0],TempB[1][1],TempB[1][2],TempB[2][0],TempB[2][1],TempB[2][2]);
		MATRIX_Inv(TempB, invB, m);
         //printf("invB:%4f,%4f,%4f\n,%4f,%4f,%4f\n,%4f,%4f,%4f\n",invB[0][0],invB[0][1],invB[0][2],invB[1][0],invB[1][1],invB[1][2],invB[2][0],invB[2][1],invB[2][2]);
		MATRIX_Mul(transA, invB ,C, n, m, m);
        
	}
	else {
		MATRIX_Mul(transA, A, TempA, n, m, n);
		MATRIX_Inv(TempA, invA, n);
		MATRIX_Mul(invA, transA, C, m, n, n);
	}
}
// 打印矩阵
void MATRIX_Display(double *A, int n,int m)
{
    int i,j;
    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            printf("%4f  ",A[i*m+j]);
        }
        printf("\n ");
    }
}
void PD_Controller(double *Theta_Err, double *dTheta_Err,double *TAU){
    double Kp_U[3],Kp_D[3];
    double Kd_U[3],Kd_D[3];
    Kp_U[0] = 300.0;
	Kp_U[1] = 200.0;
	Kp_U[2] = 150.0;
    Kp_D[0] = 30.0;
	Kp_D[1] = 200.0;
	Kp_D[2] = 150.0;
    Kd_U[0] = 20.0;
	Kd_U[1] = 15.0;
	Kd_U[2] = 10.0;
    Kd_D[0] = 20.0;
	Kd_D[1] = 15.0;
	Kd_D[2] = 10.0;
    TAU[0] = Kp_U[0]*Theta_Err[0]+Kd_U[0]*dTheta_Err[0];
    TAU[1] = Kp_U[1]*Theta_Err[1]+Kd_U[1]*dTheta_Err[1];
    TAU[2] = Kp_U[2]*Theta_Err[2]+Kd_U[2]*dTheta_Err[2];
    TAU[3] = Kp_D[0]*Theta_Err[3]+Kd_D[0]*dTheta_Err[3];
    TAU[4] = Kp_D[1]*Theta_Err[4]+Kd_D[1]*dTheta_Err[4];
    TAU[5] = Kp_D[2]*Theta_Err[5]+Kd_D[2]*dTheta_Err[5];
}
/*************************************************************************/
/*//End of math tools*/
/*************************************************************************/
/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 * See matlabroot/simulink/src/sfuntmpl_doc.c for more details.
 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
ssSetNumSampleTimes(S, 1):  Specify the number of sample times that an S-Function block has.

 */
static void mdlInitializeSizes(SimStruct *S)
{
    /* See sfuntmpl_doc.c for more details on the macros below */

    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters 输入的变量个数*/
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }
    // 设置连续状态的个数，默认值是0
    ssSetNumContStates(S, 0);
	// 设置离散状态的个数
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 19)) return;
	/* 第0个输入变量的维数为1 */
	/* Dimension of the first input variable is 0 */
    ssSetInputPortWidth(S, 0, 1);
	// 设置输入端口的信号是否mdlOutputs函数中使用，这儿设置为true。 这里也需要注意，对于每个变量都需要设置
	ssSetInputPortDirectFeedThrough(S, 0, 1);

	ssSetInputPortWidth(S, 1, 1);
	ssSetInputPortDirectFeedThrough(S, 1, 1);

	ssSetInputPortWidth(S, 2, 1);
	ssSetInputPortDirectFeedThrough(S, 2, 1); 
	ssSetInputPortWidth(S, 3, 1);
	ssSetInputPortDirectFeedThrough(S, 3, 1);
	ssSetInputPortWidth(S, 4, 1);
	ssSetInputPortDirectFeedThrough(S, 4, 1);

	ssSetInputPortWidth(S, 5, 1);
	ssSetInputPortDirectFeedThrough(S, 5, 1);

	ssSetInputPortWidth(S, 6, 3);
	ssSetInputPortDirectFeedThrough(S, 6, 1);
	ssSetInputPortWidth(S, 7, 3);
	ssSetInputPortDirectFeedThrough(S, 7, 1);

	ssSetInputPortWidth(S, 8, 3);
	ssSetInputPortDirectFeedThrough(S, 8, 1);
	ssSetInputPortWidth(S, 9, 3);
	ssSetInputPortDirectFeedThrough(S, 9, 1);

	ssSetInputPortWidth(S, 10, 3);
	ssSetInputPortDirectFeedThrough(S, 10, 1);
	ssSetInputPortWidth(S, 11, 3);
	ssSetInputPortDirectFeedThrough(S, 11, 1);
	ssSetInputPortWidth(S, 12, 9);
	ssSetInputPortDirectFeedThrough(S, 12, 1);
	ssSetInputPortWidth(S, 13, 3);
	ssSetInputPortDirectFeedThrough(S, 13, 1);
	ssSetInputPortWidth(S, 14, 4);
	ssSetInputPortDirectFeedThrough(S, 14, 1);
    ssSetInputPortWidth(S, 15, 1);
	ssSetInputPortDirectFeedThrough(S, 15, 1);
    ssSetInputPortWidth(S, 16, 5);
	ssSetInputPortDirectFeedThrough(S, 16, 1);
    ssSetInputPortWidth(S, 17, 6);
	ssSetInputPortDirectFeedThrough(S, 17, 1);
    ssSetInputPortWidth(S, 18, 13);
	ssSetInputPortDirectFeedThrough(S, 18, 1);


    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
     */
   
/*
    if (!ssSetNumOutputPorts(S, 17)) return;
    
	ssSetOutputPortWidth(S, 0, 1);
	ssSetOutputPortWidth(S, 1, 1);
	ssSetOutputPortWidth(S, 2, 1); 
	ssSetOutputPortWidth(S, 3, 1);
	ssSetOutputPortWidth(S, 4, 1);
	ssSetOutputPortWidth(S, 5, 3);
	ssSetOutputPortWidth(S, 6, 3);
	ssSetOutputPortWidth(S, 7, 3);
	ssSetOutputPortWidth(S, 8, 4); 
    ssSetOutputPortWidth(S, 9, 1);
	ssSetOutputPortWidth(S, 10, 1);
    ssSetOutputPortWidth(S, 11, 1);
    ssSetOutputPortWidth(S, 12, 1);
    ssSetOutputPortWidth(S, 13, 1);
    ssSetOutputPortWidth(S, 14, 1);
    ssSetOutputPortWidth(S, 15, 1);
    ssSetOutputPortWidth(S, 16, 13);
*/	

if (!ssSetNumOutputPorts(S, 16)) return; //设置输出变量的个数
    
	ssSetOutputPortWidth(S, 0, 1);
	ssSetOutputPortWidth(S, 1, 1);
	ssSetOutputPortWidth(S, 2, 1); 
	ssSetOutputPortWidth(S, 3, 1);
	ssSetOutputPortWidth(S, 4, 1);
	ssSetOutputPortWidth(S, 5, 3);
	ssSetOutputPortWidth(S, 6, 3);
	//ssSetOutputPortWidth(S, 7, 3);
	ssSetOutputPortWidth(S, 7, 4); 
    ssSetOutputPortWidth(S, 8, 1);
	ssSetOutputPortWidth(S, 9, 1);
    ssSetOutputPortWidth(S, 10, 1);
    ssSetOutputPortWidth(S, 11, 1);
    ssSetOutputPortWidth(S, 12, 1);
    ssSetOutputPortWidth(S, 13, 1);
    ssSetOutputPortWidth(S, 14, 1);
    ssSetOutputPortWidth(S, 15, 13);
	
    ssSetNumSampleTimes(S, 1);  //设置采样时间，此处为1s。
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Specify the sim state compliance to be same as a built-in block */
    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

    ssSetOptions(S, 0);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
//#define CONTINUOUS_SAMPLE_TIME 0.025
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, 0.005);
    ssSetOffsetTime(S, 0, 0.0);

}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {
    
    
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	// The properties of the bodies, including length, mass and inertia.
    static real_T b0_U = 0.3*1.414213526;
	static real_T b0_D = 0.3*1.414213526;
	static real_T a1_U = 0.2;
	static real_T a2_U = 0.2;
	static real_T a3_U = 0.2;
	static real_T b1_U = 0.2;
	static real_T b2_U = 0.2;
	static real_T b3_U = 0.2;
	static real_T a1_D = 0.2;
	static real_T a2_D = 0.2;
	static real_T a3_D = 0.2;
	static real_T b1_D = 0.2;
	static real_T b2_D = 0.2;
	static real_T b3_D = 0.2;
	static real_T l1_U = 0.4;
	static real_T l2_U = 0.4;
	static real_T l3_U = 0.4;
	static real_T l1_D = 0.4;
	static real_T l2_D = 0.4;
	static real_T l3_D = 0.4;

	static real_T m_B = 44.0;
	static real_T I_B = 44.0;
    
	static real_T m1_U = 2.49;
	static real_T I1_U = 0.208;

	static real_T m2_U = 2.49;
	static real_T I2_U = 0.208;

	static real_T m3_U = 2.49;
	static real_T I3_U = 0.208;


	static real_T m1_D = 2.49;
	static real_T I1_D = 0.208;

	static real_T m2_D = 2.49;
	static real_T I2_D = 0.208;

	static real_T m3_D = 2.49;
	static real_T I3_D = 0.208;

	static real_T M =  58.94;
	static real_T fai = 3.14159265/4.0;
    
    real_T Je_U_11;
	real_T Je_U_21;
	real_T Je_U_12;
	real_T Je_U_22;
	real_T Je_U_13;
	real_T Je_U_23;
	real_T Je_U_14;
	real_T Je_U_24;
	real_T Je_D_11;
	real_T Je_D_21;
	real_T Je_D_12;
	real_T Je_D_22;
	real_T Je_D_13;
	real_T Je_D_23;
	real_T Je_D_14;
	real_T Je_D_24;

	real_T H_PBthx_U;
	real_T H_PBthx_D;
	real_T H_PBthx;

	real_T H_PBthy_U;
	real_T H_PBthy_D;
	real_T H_PBthy;



	real_T H_PDM1x;
	real_T H_PDM1y;

	real_T H_PDM2x;
	real_T H_PDM2y;

	real_T H_PDM3x;
	real_T H_PDM3y;


	real_T H_PUM1x;
	real_T H_PUM1y;

	real_T H_PUM2x;
	real_T H_PUM2y;

	real_T H_PUM3x;
	real_T H_PUM3y;

	real_T H_LBth0_U;
	real_T H_LBth0_D;
	real_T H_LBth0;
	real_T H_th0;

	real_T H_LLMth1_U;
	real_T H_LLMth1_D;
	real_T H_LLMth2_U;
	real_T H_LLMth2_D;
	real_T H_LLMth3_U;
	real_T H_LLMth3_D;
	real_T L_am;//系统总的角动量
    real_T Px;//X方向线动量
    real_T Py;//Y方向线动量
    real_T L;//系统角动量
    
    real_T Px_t_1;//上一步X方向线动量
    real_T Py_t_1;//上一步Y方向线动量
    real_T L_t_1;//上一步系统角动量
    real_T dthetaE_t_1;//上一步系统角动量
    real_T PcR_t_1;//上一步系统角动量

	real_T r0_0[3][1];
	real_T r1_U[3][1];
	real_T r2_U[3][1];
	real_T r3_U[3][1];
	real_T r1_D[3][1];
	real_T r2_D[3][1];
	real_T r3_D[3][1];
	real_T r_T[3][1];
	real_T r_g[3][1];
	real_T Je_U[2][4];
	real_T Je_D[2][4];
	real_T H_PB[2][3];
	real_T H_PDM[2][3];
	real_T H_PUM[2][3];
	
	real_T re_U[3][1];
	real_T re_D[3][1];
    real_T re_the_U_t_1[3][1];
	
	real_T v0_0[3][1];
	real_T v1_U[3][1];
	real_T v2_U[3][1];
	real_T v3_U[3][1];
	real_T v1_D[3][1];
	real_T v2_D[3][1];
	real_T v3_D[3][1];
	real_T v_T[3][1];
	real_T v_g[3][1];

	real_T ve_U[3][1];
	real_T ve_D[3][1];

    real_T ve_U_t_1[3][1];

	real_T H_LB[1][3];
	real_T H_LLM_U[1][3];  //为什么要定义成这种形式？不能直接定义成 H_LLM_U[3] 吗？
	real_T H_LLM_D[1][3];
	real_T H_U[1][3];
	real_T H_D[1][3];
	real_T Je_U_ac[2][3];
	real_T Je_D_ac[2][3];
	real_T Jep_UU[6][6];
    real_T B[3][2];
    real_T A[2][3];
    real_T D[2][2];
    
    real_T cp_D[2][2];
    real_T va_D[2][2];
	//real_T H_PDM_ac[2][3];
	//real_T H_PUM_ac[2][3];
    real_T theta0;
	real_T theta0_U;
	real_T theta0_D;
	real_T theta1_U;
	real_T theta2_U;
	real_T theta3_U;
    real_T thetaE_U;
	real_T theta1_D;
	real_T theta2_D;
	real_T theta3_D;

	real_T dtheta0;
	real_T dtheta0_U;
	real_T dtheta0_D;
	real_T dtheta1_U;
	real_T dtheta2_U;
	real_T dtheta3_U;
	real_T dtheta1_D;
	real_T dtheta2_D;
	real_T dtheta3_D;
    real_T Angle_Err[6];
    real_T Angular_Err[6];
    real_T Tau_Compulating[6];
    
    real_T H_UD[1][6];
    real_T H_UD_Pinv[6][1];
    real_T C_E_identity[6][6];
    real_T C_eta[6][1];
    real_T dTheta_UD[6][1];
    real_T C_Q_t_1[1][1];
    real_T C_Q_t[1][1];
    real_T omega_0[1][1];
    real_T C_W_t[1][6];
    real_T C_W_t_1[1][6];
    real_T C_N_t[1][1];
    real_T lambda[1][1];
	real_T y_t[6][1];
	real_T Temp1[6][1];
	real_T Temp2[6][6];
	real_T Temp3[6][6];
	real_T Temp4[6][1];
	
	real_T Temp5[1][1];
	real_T Temp6[1][1];
	real_T Temp7[1][1];
	real_T Temp8[1][3];
	real_T Temp9[1][3];
	real_T Tempa[3][1];
	real_T Tempb[1][1];
	real_T Tempc[1][6];
	real_T Tempd[1][6];
	real_T Tempe[1][6];
	real_T K_t[6][1];
	
	real_T omega_UD_dst[6][1];
	real_T L_m[1][1];
	real_T C_dTheta_U_RL_t[6][1];
	real_T C_dTheta_U_RL_t_1[6][1];
    real_T Para_Iden_y[3][1];
    real_T Para_Iden_x[3][1];
    real_T Para_Iden_A[3][3];
    real_T Para_Iden_A_Inv[3][3];
    real_T dtheta_T;
    real_T dPIA_01;
    real_T dPIA_11;

    real_T Mq[9][9];
    real_T Cq[9][1];
    real_T C[9][9];
	time_T t = ssGetT(S);
	/*********** InputRealPtrsType ssGetInputPortRealSignalPtrs(SimStruct *S,int_T port) ***************/
	/***
		S： SimStruct representing an S-Function block.
		port: Index of the port whose signal is required.
	***/
	// Get pointers to signals of type double connected to an input port.
    InputRealPtrsType 	 Base_x_in = ssGetInputPortRealSignalPtrs( S, 0 );  //基座X方向位移
    InputRealPtrsType 	 Base_y_in = ssGetInputPortRealSignalPtrs( S, 1 ); //基座Y方向位移
    InputRealPtrsType 	 Theta0_in = ssGetInputPortRealSignalPtrs( S, 2 ); //基座绕Z轴转动
    InputRealPtrsType 	 Base_dx_in = ssGetInputPortRealSignalPtrs( S, 3 ); //基座X方向速度
    InputRealPtrsType 	 Base_dy_in = ssGetInputPortRealSignalPtrs( S, 4 ); //基座Y方向速度
    InputRealPtrsType 	 dTheta0_in = ssGetInputPortRealSignalPtrs( S, 5 ); //基座绕Z轴转动角速度
    InputRealPtrsType 	 Theta_U_in = ssGetInputPortRealSignalPtrs( S, 6 ); // U杆三个关节的关节角度
    InputRealPtrsType 	 dTheta_U_in = ssGetInputPortRealSignalPtrs( S, 7 ); // U杆三个关节的关节角速度
    InputRealPtrsType 	 ddTheta_U_in = ssGetInputPortRealSignalPtrs( S, 8 ); // U杆三个关节的关节角加速度
    InputRealPtrsType 	 Theta_D_in = ssGetInputPortRealSignalPtrs( S, 9 );   // D杆三个关节的关节角度
    InputRealPtrsType 	 dTheta_D_in = ssGetInputPortRealSignalPtrs( S, 10 ); // D杆三个关节的关节角速度
    InputRealPtrsType 	 ddTheta_D_in = ssGetInputPortRealSignalPtrs( S, 11 ); // D杆三个关节的关节角加速度
    InputRealPtrsType 	 Link_CG_R_U_in = ssGetInputPortRealSignalPtrs( S, 12 );
    InputRealPtrsType 	 Link_CG_dR_U_in = ssGetInputPortRealSignalPtrs( S, 13 );
    InputRealPtrsType 	 Link_CG_ddR_U_in = ssGetInputPortRealSignalPtrs( S, 14 );  
    InputRealPtrsType 	 Link_CG_R_D_in = ssGetInputPortRealSignalPtrs( S, 15 );  // Target 的角速度
    InputRealPtrsType 	 P_L_Mon_t_1 = ssGetInputPortRealSignalPtrs( S, 16 );  //Output (0-4), t-1时刻基座位姿和U臂的角度 角速度
    InputRealPtrsType 	 Link_CG_ddR_D_in = ssGetInputPortRealSignalPtrs( S, 17 ); // 末端速度和角速度
    InputRealPtrsType 	 Iteration_in = ssGetInputPortRealSignalPtrs( S, 18 );

   
    // Get a pointer to an output signal of type double (real_T)
     real_T * Base_x_out = ssGetOutputPortRealSignal( S, 0 );
     real_T * Base_y_out = ssGetOutputPortRealSignal( S, 1 );
     real_T * Theta0_out = ssGetOutputPortRealSignal( S, 2 );
     real_T * Theta_U_out = ssGetOutputPortRealSignal( S, 3 );  // 末端角速度？
     real_T * dTheta_U_out = ssGetOutputPortRealSignal( S, 4 );
     real_T * ddTheta_U_out = ssGetOutputPortRealSignal( S, 5 );
     real_T * Theta_D_out = ssGetOutputPortRealSignal( S, 6 );
    // real_T * dTheta_D_out = ssGetOutputPortRealSignal( S, 7 );
     real_T * ddTheta_D_out = ssGetOutputPortRealSignal( S, 7 );
     real_T * Tau1_U = ssGetOutputPortRealSignal( S, 8 );
     real_T * Tau2_U = ssGetOutputPortRealSignal( S, 9 );
     real_T * Tau3_U = ssGetOutputPortRealSignal( S, 10 );
     real_T * Tau1_D = ssGetOutputPortRealSignal( S, 11 );
     real_T * Tau2_D = ssGetOutputPortRealSignal( S, 12 );
     real_T * Tau3_D = ssGetOutputPortRealSignal( S, 13 );
     real_T * Counter = ssGetOutputPortRealSignal( S, 14 );
     real_T * Iteration_out = ssGetOutputPortRealSignal( S,15 );
     
     theta0 = * Theta0_in[0];

     theta0_U = theta0;
     theta0_D = theta0;
     
	theta1_U = * Theta_U_in[0];
	theta2_U = * Theta_U_in[1];
	theta3_U = * Theta_U_in[2];
	theta1_D = * Theta_D_in[0];
	theta2_D = * Theta_D_in[1];
	theta3_D = * Theta_D_in[2];

	dtheta0 = * dTheta0_in[0];

	dtheta1_U = * dTheta_U_in[0];
	dtheta2_U = * dTheta_U_in[1];
	dtheta3_U = * dTheta_U_in[2];
	dtheta1_D = * dTheta_D_in[0];
	dtheta2_D = * dTheta_D_in[1];
	dtheta3_D = * dTheta_D_in[2];
    r0_0[0][0] = * Base_x_in[0];
    r0_0[1][0] = * Base_y_in[0];
    v0_0[0][0] = * Base_dx_in[0];
    v0_0[1][0] = * Base_dy_in[0];
    C_dTheta_U_RL_t_1[0][0] = * Iteration_in[0];
    C_dTheta_U_RL_t_1[1][0] = * Iteration_in[1];
    C_dTheta_U_RL_t_1[2][0] = * Iteration_in[2];
    C_dTheta_U_RL_t_1[3][0] = * Iteration_in[3]; 
    C_dTheta_U_RL_t_1[4][0] = * Iteration_in[4]; 
    C_dTheta_U_RL_t_1[5][0] = * Iteration_in[5]; 
    C_Q_t_1[0][0] = * Iteration_in[6]; 
    C_W_t_1[0][0] = * Iteration_in[7];
    C_W_t_1[0][1] = * Iteration_in[8];
    C_W_t_1[0][2] = * Iteration_in[9];
    C_W_t_1[0][3] = * Iteration_in[10];
    C_W_t_1[0][4] = * Iteration_in[11];
    C_W_t_1[0][5] = * Iteration_in[12];
    dtheta_T = * Link_CG_R_D_in[0];
   ve_U[0][0]= * Link_CG_ddR_D_in[0];  // 末端 X 方向速度
   ve_U[1][0]= * Link_CG_ddR_D_in[1];  // 末端 Y 方向速度
   ve_U[2][0]= * Link_CG_ddR_D_in[5];  // 末端绕 Z 方向角速度
   Px_t_1 =  * P_L_Mon_t_1[0];  // t-1时刻基座X方向位置
   Py_t_1 =  * P_L_Mon_t_1[1];  // t-1时刻基座Y方向位置
   L_t_1 =  * P_L_Mon_t_1[2];   // t-1时刻基座绕Z方向旋转角度
   dthetaE_t_1 =  * P_L_Mon_t_1[3];
   PcR_t_1 =  * P_L_Mon_t_1[4];
   ve_U_t_1[0][0] =  * Link_CG_ddR_U_in[0];
   ve_U_t_1[1][0] =  * Link_CG_ddR_U_in[1];
   dPIA_01 =  * Link_CG_ddR_U_in[2];
   dPIA_11 =  * Link_CG_ddR_U_in[3];
   re_the_U_t_1[0][0] =  * Link_CG_dR_U_in[0];
   re_the_U_t_1[1][0] =  * Link_CG_dR_U_in[1];
   re_the_U_t_1[2][0] =  * Link_CG_dR_U_in[2];
   
   lambda[0][0] = 0.9999;
   
  // printf("r0:%4f,%4f\n",r0_0[0][0],r0_0[1][0]);
    r1_U[0][0] = r0_0[0][0]+b0_U*cos(theta0_U+fai)+a1_U*cos(theta0_U+theta1_U);
	r1_U[1][0] = r0_0[1][0]+b0_U*sin(theta0_U+fai)+a1_U*sin(theta0_U+theta1_U);
	//printf("r1:%4f,%4f\n",r1_U[0][0],r1_U[1][0]);
	r2_U[0][0] = r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U);
	r2_U[1][0] = r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U);
	//printf("r2:%4f,%4f\n",r2_U[0][0],r2_U[1][0]);
	r3_U[0][0] = r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U);
	r3_U[1][0] = r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U);
	//printf("r3:%4f,%4f\n",r3_U[0][0],r3_U[1][0]);
    re_U[0][0] = r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U);
	re_U[1][0] = r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U);
    //ve_U_t_1[0][0] = v0_0[0][0]-(b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta0_U-(l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta1_U-(l2_U*sin(theta0_U+theta1_U+theta2_U)+l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta2_U-(l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta3_U;
    //ve_U_t_1[1][0] = v0_0[1][0]+(b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta0_U+(l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta1_U+(l2_U*cos(theta0_U+theta1_U+theta2_U)+l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta2_U+(l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*dtheta3_U;
	r1_D[0][0] = r0_0[0][0]+b0_D*cos(theta0_D-fai)+a1_D*cos(theta0_D+theta1_D);
	r1_D[1][0] = r0_0[1][0]+b0_D*sin(theta0_D-fai)+a1_D*sin(theta0_D+theta1_D);
	//printf("r1d:%4f,%4f\n",r1_D[0][0],r1_D[1][0]);
	r2_D[0][0] = r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D);
	r2_D[1][0] = r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D);
	//printf("r2d:%4f,%4f\n",r2_D[0][0],r2_D[1][0]);
	r3_D[0][0] = r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D);
	r3_D[1][0] = r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D);
	//printf("r3d:%4f,%4f\n",r3_D[0][0],r3_D[1][0]);
    r_g[0][0] = (m_B*r0_0[0][0]+m1_U*r1_U[0][0]+m2_U*r2_U[0][0]+m3_U*r3_U[0][0]+m1_D*r1_D[0][0]+m2_D*r2_D[0][0]+m3_D*r3_D[0][0])/M;
    r_g[1][0] = (m_B*r0_0[1][0]+m1_U*r1_U[1][0]+m2_U*r2_U[1][0]+m3_U*r3_U[1][0]+m1_D*r1_D[1][0]+m2_D*r2_D[1][0]+m3_D*r3_D[1][0])/M;
   // printf("U:%4f,%4f,%4f\n",theta1_U*180.0/3.1415926,theta2_U*180.0/3.1415926,theta3_U*180.0/3.1415926);  
   // printf("U:%4f,%4f,%4f\n",theta1_D*180.0/3.1415926,theta2_D*180.0/3.1415926,theta3_D*180.0/3.1415926);  
   /***************************** 关节Jacobian矩阵 **********************************************/
   /* Je为关节Jacobian矩阵，其中Je_U_12，Je_U_22，Je_D_12，Je_D_22 几项是不是应该去掉含l1_U的项？Je0为底座Jacobian矩阵  */
    Je_U_11 = -1.0*(b0_U*sin(theta0_U+fai) + l1_U*sin(theta0_U+theta1_U) + l2_U*sin(theta0_U+theta1_U+theta2_U) + l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_21 =      (b0_U*cos(theta0_U+fai) + l1_U*cos(theta0_U+theta1_U) + l2_U*cos(theta0_U+theta1_U+theta2_U) + l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_12 = -1.0*(l1_U*sin(theta0_U+theta1_U) + l2_U*sin(theta0_U+theta1_U+theta2_U) + l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_22 =      (l1_U*cos(theta0_U+theta1_U) + l2_U*cos(theta0_U+theta1_U+theta2_U) + l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_13 = -1.0*(l2_U*sin(theta0_U+theta1_U+theta2_U) + l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_23 =      (l2_U*cos(theta0_U+theta1_U+theta2_U) + l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_14 = -1.0*(l3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U));
	Je_U_24 =      (l3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U));

	
	Je_D_11 = -1.0*(b0_D*sin(theta0_D-fai) + l1_D*sin(theta0_D+theta1_D) + l2_D*sin(theta0_D+theta1_D+theta2_D) + l3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_21 =      (b0_D*cos(theta0_D-fai) + l1_D*cos(theta0_D+theta1_D) + l2_D*cos(theta0_D+theta1_D+theta2_D) + l3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_12 = -1.0*(l1_D*sin(theta0_D+theta1_D) + l2_D*sin(theta0_D+theta1_D+theta2_D) + l3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_22 =      (l1_D*cos(theta0_D+theta1_D) + l2_D*cos(theta0_D+theta1_D+theta2_D) + l3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_13 = -1.0*(l2_D*sin(theta0_D+theta1_D+theta2_D) + l3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_23 =      (l2_D*cos(theta0_D+theta1_D+theta2_D) + l3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_14 = -1.0*(l3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D));
	Je_D_24 =      (l3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D));
	
	/***************************** 线动量守恒 **********************************************/
	// X方向线动量基座项系数
	H_PBthx_U = -1.0*(m1_U*(b0_U*sin(theta0_U+fai)+a1_U*sin(theta0_U+theta1_U)) + m2_U*(b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U)) + m3_U*(b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)));
	H_PBthx_D = -1.0*(m1_D*(b0_D*sin(theta0_D-fai)+a1_D*sin(theta0_D+theta1_D)) + m2_D*(b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D)) + m3_D*(b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)));
	H_PBthx = H_PBthx_U + H_PBthx_D;

	//printf("H_PBthx:%4f\n",H_PBthx);
	// Y方向线动量基座项系数
	H_PBthy_U = (m1_U*(b0_U*cos(theta0_U+fai)+a1_U*cos(theta0_U+theta1_U)) + m2_U*(b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U)) + m3_U*(b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U)));
	H_PBthy_D = (m1_D*(b0_D*cos(theta0_D-fai)+a1_D*cos(theta0_D+theta1_D)) + m2_D*(b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D)) + m3_D*(b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D)));
	H_PBthy = H_PBthy_U + H_PBthy_D;

	//printf("H_PBthy:%4f\n",H_PBthy);
    // U臂线动量关节系数
	H_PUM1x = -1.0*(m1_U*(a1_U*sin(theta0_U+theta1_U)) + m2_U*(l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U)) + m3_U*(l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)));
	H_PUM1y =      (m1_U*(a1_U*cos(theta0_U+theta1_U)) + m2_U*(l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U)) + m3_U*(l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U)));

	H_PUM2x = -1.0*(m2_U*(a2_U*sin(theta0_U+theta1_U+theta2_U)) + m3_U*(l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)));
	H_PUM2y =      (m2_U*(a2_U*cos(theta0_U+theta1_U+theta2_U)) + m3_U*(l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U)));

	H_PUM3x = -1.0*(m3_U*(a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)));
	H_PUM3y =      (m3_U*(a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U)));
	//printf("H_PUM1x:%4f\n",H_PUM1x);
	//printf("H_PUM1y:%4f\n",H_PUM1y);
	//printf("H_PUM2x:%4f\n",H_PUM2x);
	//printf("H_PUM2y:%4f\n",H_PUM2y);
	//printf("H_PUM3x:%4f\n",H_PUM3x);
	//printf("H_PUM3y:%4f\n",H_PUM3y);

	// D臂线动量关节系数
	H_PDM1x = -1.0*(m1_D*(a1_D*sin(theta0_D+theta1_D)) + m2_D*(l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D)) + m3_D*(l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)));
	H_PDM1y =      (m1_D*(a1_D*cos(theta0_D+theta1_D)) + m2_D*(l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D)) + m3_D*(l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D)));

	H_PDM2x = -1.0*(m2_D*(a2_D*sin(theta0_D+theta1_D+theta2_D)) + m3_D*(l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)));
	H_PDM2y =      (m2_D*(a2_D*cos(theta0_D+theta1_D+theta2_D)) + m3_D*(l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D)));

	H_PDM3x = -1.0*(m3_D*(a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)));
	H_PDM3y = (m3_D*(a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D)));

	//printf("H_PDM1x:%4f\n",H_PDM1x);
	//printf("H_PDM1y:%4f\n",H_PDM1y);
	//printf("H_PDM2x:%4f\n",H_PDM2x);
	//printf("H_PDM2y:%4f\n",H_PDM2y);
	//printf("H_PDM3x:%4f\n",H_PDM3x);
	//printf("H_PDM3y:%4f\n",H_PDM3y);
	
	// X Y方向上线动量
    Px = M*v0_0[0][0] + H_PBthx*dtheta0 + H_PUM1x*dtheta1_U + H_PUM2x*dtheta2_U + H_PUM3x*dtheta3_U + H_PDM1x*dtheta1_D + H_PDM2x*dtheta2_D + H_PDM3x*dtheta3_D;
	Py = M*v0_0[1][0] + H_PBthy*dtheta0 + H_PUM1y*dtheta1_U + H_PUM2y*dtheta2_U + H_PUM3y*dtheta3_U + H_PDM1y*dtheta1_D + H_PDM2y*dtheta2_D + H_PDM3y*dtheta3_D;

	/***************************** 角动量守恒 **********************************************/
	// 基座角动量系数
	H_LBth0_U  = I_B + I1_U + I2_U + I3_U + m3_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*(b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*(b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)))+ 
											m2_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U))*(b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U))*(b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U)))+
											m1_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+a1_U*cos(theta0_U+theta1_U))*(b0_U*cos(theta0_U+fai)+a1_U*cos(theta0_U+theta1_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+a1_U*sin(theta0_U+theta1_U))*(b0_U*sin(theta0_U+fai)+a1_U*sin(theta0_U+theta1_U)));

	H_LBth0_D =        I1_D + I2_D + I3_D + m3_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))*(b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D))*(b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)))+
											m2_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D))*(b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D))*(b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D)))+ 
											m1_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+a1_D*cos(theta0_D+theta1_D))*(b0_D*cos(theta0_D-fai)+a1_D*cos(theta0_D+theta1_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+a1_D*sin(theta0_D+theta1_D))*(b0_D*sin(theta0_D-fai)+a1_D*sin(theta0_D+theta1_D)));

	H_LBth0 =  H_LBth0_U + H_LBth0_D;
	// 双臂角动量系数
	H_LLMth1_U =  I1_U +  I2_U + I3_U + m3_U*((r0_0[0][0] + b0_U*cos(theta0_U+fai) + l1_U*cos(theta0_U+theta1_U) + l2_U*cos(theta0_U+theta1_U+theta2_U) + a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*(l1_U*cos(theta0_U+theta1_U) + l2_U*cos(theta0_U+theta1_U+theta2_U) + a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U)) + (r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U) + l2_U*sin(theta0_U+theta1_U+theta2_U) + a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*(l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)))+
										m2_U*((r0_0[0][0] + b0_U*cos(theta0_U+fai) + l1_U*cos(theta0_U+theta1_U) + a2_U*cos(theta0_U+theta1_U+theta2_U))*(l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U)) + (r0_0[1][0] + b0_U*sin(theta0_U+fai) + l1_U*sin(theta0_U+theta1_U) + a2_U*sin(theta0_U+theta1_U+theta2_U))*(l1_U*sin(theta0_U+theta1_U) + a2_U*sin(theta0_U+theta1_U+theta2_U)))+
										m1_U*((r0_0[0][0] + b0_U*cos(theta0_U+fai) + a1_U*cos(theta0_U+theta1_U))*(a1_U*cos(theta0_U+theta1_U)) + (r0_0[1][0]+b0_U*sin(theta0_U+fai)+a1_U*sin(theta0_U+theta1_U))*(a1_U*sin(theta0_U+theta1_U)));

	H_LLMth1_D = I1_D + I2_D + I3_D  +  m3_D*((r0_0[0][0] + b0_D*cos(theta0_D-fai) + l1_D*cos(theta0_D+theta1_D) + l2_D*cos(theta0_D+theta1_D+theta2_D) + a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))*(l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D) + a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D))*(l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)))+
										m2_D*((r0_0[0][0] + b0_D*cos(theta0_D-fai) + l1_D*cos(theta0_D+theta1_D) + a2_D*cos(theta0_D+theta1_D+theta2_D))*(l1_D*cos(theta0_D+theta1_D) + a2_D*cos(theta0_D+theta1_D+theta2_D)) + (r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D))*(l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D)))+
										m1_D*((r0_0[0][0] + b0_D*cos(theta0_D-fai) + a1_D*cos(theta0_D+theta1_D))*(a1_D*cos(theta0_D+theta1_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+a1_D*sin(theta0_D+theta1_D))*(a1_D*sin(theta0_D+theta1_D)));

	H_LLMth2_U = I2_U + I3_U  + m3_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*(l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*(l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)))+ 
								m2_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+a2_U*cos(theta0_U+theta1_U+theta2_U))*(a2_U*cos(theta0_U+theta1_U+theta2_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+a2_U*sin(theta0_U+theta1_U+theta2_U))*(a2_U*sin(theta0_U+theta1_U+theta2_U)));

	H_LLMth2_D = I2_D + I3_D  + m3_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))*(l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D))*(l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)))+
								m2_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+a2_D*cos(theta0_D+theta1_D+theta2_D))*(a2_D*cos(theta0_D+theta1_D+theta2_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+a2_D*sin(theta0_D+theta1_D+theta2_D))*(a2_D*sin(theta0_D+theta1_D+theta2_D)));

	H_LLMth3_U = I3_U + m3_U*((r0_0[0][0]+b0_U*cos(theta0_U+fai)+l1_U*cos(theta0_U+theta1_U)+l2_U*cos(theta0_U+theta1_U+theta2_U)+a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))*(a3_U*cos(theta0_U+theta1_U+theta2_U+theta3_U))+(r0_0[1][0]+b0_U*sin(theta0_U+fai)+l1_U*sin(theta0_U+theta1_U)+l2_U*sin(theta0_U+theta1_U+theta2_U)+a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U))*(a3_U*sin(theta0_U+theta1_U+theta2_U+theta3_U)));

	H_LLMth3_D = I3_D + m3_D*((r0_0[0][0]+b0_D*cos(theta0_D-fai)+l1_D*cos(theta0_D+theta1_D)+l2_D*cos(theta0_D+theta1_D+theta2_D)+a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))*(a3_D*cos(theta0_D+theta1_D+theta2_D+theta3_D))+(r0_0[1][0]+b0_D*sin(theta0_D-fai)+l1_D*sin(theta0_D+theta1_D)+l2_D*sin(theta0_D+theta1_D+theta2_D)+a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D))*(a3_D*sin(theta0_D+theta1_D+theta2_D+theta3_D)));

    // 基座线动量系数矩阵
    H_PB[0][0] = M; 
	H_PB[0][1] = 0.0; 
	H_PB[0][2] = H_PBthx;
    H_PB[1][0] = 0.0; 
	H_PB[1][1] = M; 
	H_PB[1][2] = H_PBthy;
	// D臂线动量系数矩阵
    H_PDM[0][0] = H_PDM1x; 
	H_PDM[0][1] = H_PDM2x; 
	H_PDM[0][2] = H_PDM3x;
    H_PDM[1][0] = H_PDM1y; 
	H_PDM[1][1] = H_PDM2y; 
	H_PDM[1][2] = H_PDM3y;
	// U臂线动量系数矩阵
    H_PUM[0][0] = H_PUM1x;
	H_PUM[0][1] = H_PUM2x;
	H_PUM[0][2] = H_PUM3x;
    H_PUM[1][0] = H_PUM1y;
	H_PUM[1][1] = H_PUM2y;
	H_PUM[1][2] = H_PUM3y;
/*   
    H_PDM_ac[0][0] = H_PDM1x/M;
	H_PDM_ac[0][1] = H_PDM2x/M;
	H_PDM_ac[0][2] = H_PDM3x/M;
    H_PDM_ac[1][0] = H_PDM1y/M;
	H_PDM_ac[1][1] = H_PDM2y/M;
	H_PDM_ac[1][2] = H_PDM3y/M;
	
    H_PUM_ac[0][0] = H_PUM1x/M;
	H_PUM_ac[0][1] = H_PUM2x/M;
	H_PUM_ac[0][2] = H_PUM3x/M;
    H_PUM_ac[1][0] = H_PUM1y/M;
	H_PUM_ac[1][1] = H_PUM2y/M;
	H_PUM_ac[1][2] = H_PUM3y/M;
*/
	// 基座角动量系数矩阵
    H_LB[0][0] = -1.0*M*r_g[1][0];
	H_LB[0][1] = M*r_g[0][0];
	H_LB[0][2] = H_LBth0;
	// 双臂角动量系数矩阵
    H_LLM_U[0][0] = H_LLMth1_U;
	H_LLM_U[0][1] = H_LLMth2_U;
	H_LLM_U[0][2] = H_LLMth3_U;
    H_LLM_D[0][0] = H_LLMth1_D;
	H_LLM_D[0][1] = H_LLMth2_D;
	H_LLM_D[0][2] = H_LLMth3_D;
    
	// 这一部分怎么来的？？？？？？？？？？？？？？？？
//将线动量方程与角动量方程合并后，可表示为H_U*dTheta_U + H_D*dTheta_L + H_th0*dTheta0 = 0 
//H_U、H_D为1X3矩阵
    H_U[0][0] = H_LLM_U[0][0] + H_PUM[0][0]*r_g[1][0] - H_PUM[1][0]*r_g[0][0];
    H_U[0][1] = H_LLM_U[0][1] + H_PUM[0][1]*r_g[1][0] - H_PUM[1][1]*r_g[0][0];
    H_U[0][2] = H_LLM_U[0][2] + H_PUM[0][2]*r_g[1][0] - H_PUM[1][2]*r_g[0][0];
    
	H_D[0][0] = H_LLM_D[0][0] + H_PDM[0][0]*r_g[1][0] - H_PDM[1][0]*r_g[0][0];
    H_D[0][1] = H_LLM_D[0][1] + H_PDM[0][1]*r_g[1][0] - H_PDM[1][1]*r_g[0][0];
    H_D[0][2] = H_LLM_D[0][2] + H_PDM[0][2]*r_g[1][0] - H_PDM[1][2]*r_g[0][0];
    
	H_th0 = H_LBth0-(H_PBthy*r_g[0][0]-H_PBthx*r_g[1][0]);
    
    // 系统总体角动量？？？？？	
    L_am = H_th0*dtheta0 + H_U[0][0]*dtheta1_U + H_U[0][1]*dtheta2_U + H_U[0][2]*dtheta3_U + H_D[0][0]*dtheta1_D + H_D[0][1]*dtheta2_D+H_D[0][2]*dtheta3_D;
    //L = H_LB[0][0]*v0_0[0][0]+H_LB[0][1]*v0_0[1][0]+ H_LBth0*dtheta0+H_LLM_U[0][0]*dtheta1_U+H_LLM_U[0][1]*dtheta2_U+H_LLM_U[0][2]*dtheta3_U+H_LLM_D[0][0]*dtheta1_D+H_LLM_D[0][1]*dtheta2_D+H_LLM_D[0][2]*dtheta3_D + (Px)*re_U[1][0]-(Py)*re_U[0][0];
    // 线动量
	L = H_LB[0][0]*v0_0[0][0] + H_LB[0][1]*v0_0[1][0] + H_LBth0*dtheta0 + H_LLM_U[0][0]*dtheta1_U + H_LLM_U[0][1]*dtheta2_U + H_LLM_U[0][2]*dtheta3_U + H_LLM_D[0][0]*dtheta1_D + H_LLM_D[0][1]*dtheta2_D + H_LLM_D[0][2]*dtheta3_D;
    printf("L:%4f\n",H_LB[0][0]*v0_0[0][0]+H_LB[0][1]*v0_0[1][0]+ H_LBth0*dtheta0+H_LLM_U[0][0]*dtheta1_U+H_LLM_U[0][1]*dtheta2_U+H_LLM_U[0][2]*dtheta3_U+H_LLM_D[0][0]*dtheta1_D+H_LLM_D[0][1]*dtheta2_D+H_LLM_D[0][2]*dtheta3_D);
  
    MATRIX_SetZero(y_t, 6,1);
    MATRIX_SetZero(Temp1, 6,1);
    MATRIX_SetZero(Temp2, 6,6);
    MATRIX_SetZero(Temp3, 6,6);
    MATRIX_SetZero(Temp4, 6,1);
    MATRIX_SetZero(Temp5, 1,1);
    MATRIX_SetZero(Temp6, 1,1);
    MATRIX_SetZero(Temp7, 1,1);
    MATRIX_SetZero(Temp8, 1,3);
    MATRIX_SetZero(Temp9, 1,3);
    MATRIX_SetZero(Tempa, 3,1);
    MATRIX_SetZero(Tempb, 1,1);
    MATRIX_SetZero(Tempc, 1,6);
    MATRIX_SetZero(Tempd, 1,6);
    MATRIX_SetZero(Tempe, 1,6);
    MATRIX_SetZero(K_t, 6,1);
    MATRIX_SetZero(omega_0, 1,1);
    MATRIX_SetZero(omega_UD_dst, 6,1);
    printf("C_dTheta_U_RL_t_1:%4f,%4f,%4f,%4f,%4f,%4f\n",C_dTheta_U_RL_t_1[0][0],C_dTheta_U_RL_t_1[1][0],C_dTheta_U_RL_t_1[2][0],C_dTheta_U_RL_t_1[3][0],C_dTheta_U_RL_t_1[4][0],C_dTheta_U_RL_t_1[5][0]);
    printf("C_Q_t_1:%4f\n",C_Q_t_1[0][0]);
    printf("C_W_t_1:%4f,%4f,%4f,%4f,%4f,%4f\n",C_W_t_1[0][0],C_W_t_1[0][1],C_W_t_1[0][2],C_W_t_1[0][3],C_W_t_1[0][4],C_W_t_1[0][5]);
    // 关节角速度
	dTheta_UD[0][0] = dtheta1_U;
	dTheta_UD[1][0] = dtheta2_U;
	dTheta_UD[2][0] = dtheta3_U;
	dTheta_UD[3][0] = dtheta1_D;
	dTheta_UD[4][0] = dtheta2_D;
	dTheta_UD[5][0] = dtheta3_D;
	omega_0[0][0] = dtheta0;
    H_UD[0][0] = H_U[0][0];
    H_UD[0][1] = H_U[0][1];
    H_UD[0][2] = H_U[0][2];
    H_UD[0][3] = H_D[0][0];
    H_UD[0][4] = H_D[0][1];
    H_UD[0][5] = H_D[0][2];
    MATRIX_SetZero(H_UD_Pinv,6,1);
    MATRIX_SetUnit(C_E_identity,6);
        //未考虑关节角度限制
    C_eta[0][0] = 0.41;
	C_eta[1][0] = 0.41;
	C_eta[2][0] = 0.41;
    C_eta[3][0] = 0.41;
	C_eta[4][0] = 0.41;
	C_eta[5][0] = 0.41;
        //考虑关节角度限制
        // C_eta[0][0] = -0.3*2.0*theta1_U/(100.0*100.0-theta1_U*theta1_U)/(100.0*100.0-theta1_U*theta1_U);C_eta[1][0] = -0.3*2.0*theta2_U/(90.0*90.0-theta2_U*theta2_U)/(90.0*90.0-theta2_U*theta2_U);C_eta[2][0] = -0.3*2.0*theta3_U/(60.0*60.0-theta3_U*theta3_U)/(60.0*60.0-theta3_U*theta3_U);
        // C_eta[3][0] = -0.3*2.0*theta1_D/(100.0*100.0-theta1_D*theta1_D)/(100.0*100.0-theta1_D*theta1_D);C_eta[4][0] = -0.3*2.0*theta2_D/(90.0*90.0-theta2_D*theta2_D)/(90.0*90.0-theta2_D*theta2_D);C_eta[5][0] = -0.3*2.0*theta3_D/(60.0*60.0-theta3_D*theta3_D)/(60.0*60.0-theta3_D*theta3_D);
        
		/*************** 计算反作用零空间 ****************************/
		// if中为抓捕前的信息初始化
		// L_am 为系统总体角动量
		// H_UD 为 (16)-(17) 中的H_wfai矩阵
        if( (t/0.005) < 1.0){
            L_m[0][0] = L_am;
            MATRIX_Pinv(H_UD,H_UD_Pinv,1,6);
            MATRIX_Mul(H_UD_Pinv,L_m,Temp1,6,1,1);  // (19) 式中第一项
            MATRIX_Mul(H_UD_Pinv,H_UD,Temp2,6,1,6); 
            MATRIX_Sub(C_E_identity,Temp2,Temp3,6,6); // 和上式构成（17），反作用零空间
            MATRIX_Mul(Temp3,C_eta,Temp4,6,6,1); //构成（19）第二式
            MATRIX_Add(Temp1,Temp4,C_dTheta_U_RL_t_1,6,1); // 相加得到（19）式
            //计算C_W_t
            Temp5[0][0] = H_th0;
            MATRIX_Mul(H_UD_Pinv,Temp5,Temp4,6,1,1);
            MATRIX_Tran(Temp4,C_W_t_1,6,1);
            C_Q_t_1[0][0] = 500.0;  // 为什么是500？ 
        }
        else {

          /*  L_m[0][0] = L_am;
          
            MATRIX_SetZero(H_UD_Pinv,6,1);
            MATRIX_SetUnit(C_E_identity,6);

            MATRIX_Pinv(H_UD,H_UD_Pinv,1,6);
            MATRIX_Mul(H_UD_Pinv,L_m,Temp1,6,1,1);
            MATRIX_Mul(H_UD_Pinv,H_UD,Temp2,6,1,6);
            MATRIX_Sub(C_E_identity,Temp2,Temp3,6,6);
            MATRIX_Mul(Temp3,C_eta,Temp4,6,6,1);
            MATRIX_Add(Temp1,Temp4,C_dTheta_U_RL_t,6,1);*/
         }
	//计算y(t)
	MATRIX_Sub(C_dTheta_U_RL_t_1,dTheta_UD,y_t,6,1);
	//gsl_matrix_tran(y_t,Temp4);//计算y^T(t)
	//更新W_hat(t)
		/*计算N(t)*/
	MATRIX_Mul(C_Q_t_1,omega_0,Temp5,1,1,1); // （25）中，N（t）分子
	MATRIX_Tran(omega_0,Tempb,1,1);
	MATRIX_Mul(Tempb,Temp5,Temp6,1,1,1);
	MATRIX_Add(lambda,Temp6,Tempb,1,1,1);
	MATRIX_Inv(Tempb,Temp6,1);  // （25）中，N（t）分母
	MATRIX_Mul(Temp5,Temp6,C_N_t,1,1,1); // 得到N（t）
    //printf("C_N_t:%4f\n",C_N_t[0][0]);
	/*计算Q(t)*/
	//gsl_matrix_tran(Temp5,Tempb);
	//gsl_matrix_mul(C_N_t,Tempb,Temp6);
	MATRIX_Tran(omega_0,Tempb,1,1);
	MATRIX_Mul(C_N_t,Tempb,Temp6,1,1,1);
	MATRIX_Mul(Temp6,C_Q_t_1,Temp5,1,1,1);  // 式（25）中Q(t) 减数
	MATRIX_Sub(C_Q_t_1,Temp5,Temp6,1,1,1);
	MATRIX_Inv(lambda,Tempb,1);
	MATRIX_Mul(Temp6,Tempb,C_Q_t,1,1,1);  // 得到Q(t)
	/*计算W_hat(t)*/
	MATRIX_Tran(omega_0,Tempb,1,1);
	MATRIX_Mul(Tempb,C_W_t_1,Tempc,1,1,6);
	MATRIX_Tran(y_t,Tempe,6,1);
	MATRIX_Sub(Tempe,Tempc,Tempd,1,6);
	MATRIX_Mul(C_N_t,Tempd,Tempe,1,1,6);
	MATRIX_Add(C_W_t_1,Tempe,C_W_t,1,6);  // 得到 W_hat(t) */
	//更新omega_UD_dst
	MATRIX_Tran(C_W_t,K_t,1,6);
	MATRIX_Mul(K_t,omega_0,Temp1,6,1,1);
	MATRIX_Add(Temp1,dTheta_UD,omega_UD_dst,6,1);  // 式(26)
   
   /*  omega_UD_dst 规划6个关节角的角速度                   */
    C_dTheta_U_RL_t[0][0] = omega_UD_dst[0][0];
    C_dTheta_U_RL_t[1][0] = omega_UD_dst[1][0];
    C_dTheta_U_RL_t[2][0] = omega_UD_dst[2][0];
    C_dTheta_U_RL_t[3][0] = omega_UD_dst[3][0];
    C_dTheta_U_RL_t[4][0] = omega_UD_dst[4][0];
    C_dTheta_U_RL_t[5][0] = omega_UD_dst[5][0];
	
	// 规划角速度与实际角速度之差
	Angular_Err[0] =  omega_UD_dst[0][0] - dtheta1_U;
    Angular_Err[1] =  omega_UD_dst[1][0] - dtheta2_U;
    Angular_Err[2] =  omega_UD_dst[2][0] - dtheta3_U;
    
    Angular_Err[3] =  omega_UD_dst[3][0] - dtheta1_D;
    Angular_Err[4] =  omega_UD_dst[4][0] - dtheta2_D;
    Angular_Err[5] =  omega_UD_dst[5][0] - dtheta3_D;
	
	// 误差为什么要乘以0.025？  转化成角度之差！
	Angle_Err[0] = Angular_Err[0]*(0.025);
    Angle_Err[1] = Angular_Err[1]*(0.025);
    Angle_Err[2] = Angular_Err[2]*(0.025);
    
    Angle_Err[3] = Angular_Err[3]*(0.025);
    Angle_Err[4] = Angular_Err[4]*(0.025);
    Angle_Err[5] = Angular_Err[5]*(0.025);
    
    
	//C_W_t_1[0][0] = C_W_t[0][0];
	//C_Q_t_1[0][0] = C_Q_t[0][0];
	if((C_Q_t[0][0]>600.0)||(C_Q_t[0][0]<1))
    {
			C_Q_t[0][0]=500.0;
		}
    
    
   // printf("rg:%4f,%4f\n",r_g[0][0],r_g[1][0]);||(C_Q_t[0][0]<1)
    //printf("T:%4f,%4f\n",t,t/0.025);
    //Angle_Err[0] = (theta1_U - 16.0*3.14159265/180.0)*(-1.0);
   // Angle_Err[1] = (theta2_U - (-16.0)*3.14159265/180.0)*(-1.0);
   // Angle_Err[2] = (theta3_U - (-6.0)*3.14159265/180.0)*(-1.0);
    
   // Angle_Err[3] = (theta1_D - (-16.0)*3.14159265/180.0)*(-1.0);
   // Angle_Err[4] = (theta2_D - 16.0*3.14159265/180.0)*(-1.0);
   // Angle_Err[5] = (theta3_D - 6.0*3.14159265/180.0)*(-1.0);
    
   // Angular_Err[0] = (dtheta1_U - 0.0*3.14159265/180.0)*(-1.0);
   // Angular_Err[1] = (dtheta2_U - 0.0*3.14159265/180.0)*(-1.0);
   // Angular_Err[2] = (dtheta3_U - 0.0*3.14159265/180.0)*(-1.0);
    
   // Angular_Err[3] = (dtheta1_D - 0.0*3.14159265/180.0)*(-1.0);
   // Angular_Err[4] = (dtheta2_D - 0.0*3.14159265/180.0)*(-1.0);
   // Angular_Err[5] = (dtheta3_D - 0.0*3.14159265/180.0)*(-1.0);
   // if( (t/0.025) < 1.0){
        //Tau1_U[0] = 0.0;
       // Tau2_U[0] = 0.0;
        //Tau3_U[0] = 0.0;
       // Tau1_D[0] = 0.0;
        //Tau2_D[0] = 0.0;
       // Tau3_D[0] = 0.0;
    
   // }
   // else
    //{
        PD_Controller(Angle_Err, Angular_Err,Tau_Compulating);
    
        Tau1_U[0] = Tau_Compulating[0];
        Tau2_U[0] = Tau_Compulating[1];
        Tau3_U[0] = Tau_Compulating[2];
        Tau1_D[0] = Tau_Compulating[3];
        Tau2_D[0] = Tau_Compulating[4];
        Tau3_D[0] = Tau_Compulating[5];
   // }
    printf("Tau:%4f,%4f,%4f,%4f,%4f,%4f\n",Tau1_U[0],Tau2_U[0],Tau3_U[0],Tau1_D[0],Tau2_D[0],Tau3_D[0]);
    //参数辨识过程
   Para_Iden_y[0][0]= 1.0*(ve_U_t_1[0][0] - ve_U[0][0]);
   Para_Iden_y[1][0]= 1.0*(ve_U_t_1[1][0] - ve_U[1][0]);
   Para_Iden_y[2][0]= L_t_1 - L-((Px- Px_t_1)*re_U[1][0]-(Py- Py_t_1)*re_U[0][0]) ;
   thetaE_U = theta0_U + theta1_U + theta2_U + theta3_U;
   Para_Iden_A[0][0] = Px - Px_t_1;
   Para_Iden_A[0][1] = -1.0*sin(thetaE_U)*ve_U[2][0] - dPIA_01;
   Para_Iden_A[0][2] = 0.0;
   Para_Iden_A[1][0] = Py - Py_t_1;
   Para_Iden_A[1][1] = cos(thetaE_U)*ve_U[2][0] - dPIA_11;
   Para_Iden_A[1][2] = 0.0;
   Para_Iden_A[2][0] = 0.0;               
   Para_Iden_A[2][1] = (Px- Px_t_1)*(sin(thetaE_U))-(Py- Py_t_1)*(cos(thetaE_U))-  PcR_t_1*0;     
   Para_Iden_A[2][2] = (ve_U[2][0]  - dthetaE_t_1*1);//(dtheta0 + dtheta1_U + dtheta2_U + dtheta3_U - dthetaE_t_1)*0.025;
   //printf("VE:%4f,%4f,%4f\n",ve_U[2][0] ,dthetaE_t_1,Para_Iden_A[2][2]);
  // MATRIX_Display(Para_Iden_A, 3,3);
   MATRIX_Inv(Para_Iden_A,Para_Iden_A_Inv,3);
   //MATRIX_Display(Para_Iden_A_Inv, 3,3);
   MATRIX_Mul(Para_Iden_A_Inv,Para_Iden_y,Para_Iden_x,3,3,1);
   printf("Para:%4f,%4f,%4f\n",1/Para_Iden_x[0][0],Para_Iden_x[1][0],Para_Iden_x[2][0]);
   //参数辨识结束
    Base_x_out[0] =  Px;
    Base_y_out[0] =  Py;
    Theta0_out[0] = L;
    Theta_U_out[0] = ve_U[2][0] ;//dtheta0 + dtheta1_U + dtheta2_U + dtheta3_U;
    dTheta_U_out[0] = (Px)*(sin(thetaE_U))-(Py)*(cos(thetaE_U));
    ddTheta_D_out[0] =  ve_U[0][0];
    ddTheta_D_out[1] =  ve_U[1][0];
    ddTheta_D_out[2] =  -1.0*sin(thetaE_U)*ve_U[2][0];
    ddTheta_D_out[3] =  cos(thetaE_U)*ve_U[2][0];
    Iteration_out[0] = C_dTheta_U_RL_t[0][0];
    Iteration_out[1] = C_dTheta_U_RL_t[1][0];
    Iteration_out[2] = C_dTheta_U_RL_t[2][0];
    Iteration_out[3] = C_dTheta_U_RL_t[3][0];
    Iteration_out[4] = C_dTheta_U_RL_t[4][0];
    Iteration_out[5] = C_dTheta_U_RL_t[5][0];
    Iteration_out[6] = C_Q_t[0][0];
    Iteration_out[7] = C_W_t[0][0];
    Iteration_out[8] = C_W_t[0][1];
    Iteration_out[9] = C_W_t[0][2];
    Iteration_out[10] = C_W_t[0][3];
    Iteration_out[11] = C_W_t[0][4];
    Iteration_out[12] = C_W_t[0][5];
    Theta_D_out[0] = 1.0/Para_Iden_x[0][0];
    Theta_D_out[1] = Para_Iden_x[1][0];
    Theta_D_out[2] = Para_Iden_x[2][0];
    ddTheta_U_out[0] = re_U[0][0];
    ddTheta_U_out[1] = re_U[1][0];
    ddTheta_U_out[2] = thetaE_U;
    //printf("C_dTheta_U_RL_t:%4f,%4f,%4f,%4f,%4f,%4f\n",C_dTheta_U_RL_t[0][0],C_dTheta_U_RL_t[1][0],C_dTheta_U_RL_t[2][0],C_dTheta_U_RL_t[3][0],C_dTheta_U_RL_t[4][0],C_dTheta_U_RL_t[5][0]);
    //printf("C_Q_t:%4f\n",C_Q_t_1[0][0]);
   // printf("C_W_t:%4f,%4f,%4f,%4f,%4f,%4f\n",C_W_t[0][0],C_W_t[0][1],C_W_t[0][2],C_W_t[0][3],C_W_t[0][4],C_W_t[0][5]);
   printf("omega_UD_dst:%4f,%4f,%4f,%4f,%4f,%4f\n",omega_UD_dst[0][0],omega_UD_dst[1][0],omega_UD_dst[2][0],omega_UD_dst[3][0],omega_UD_dst[4][0],omega_UD_dst[5][0]);
    printf("Angle_Err:%4f,%4f,%4f,%4f,%4f,%4f\n",Angle_Err[0],Angle_Err[1],Angle_Err[2],Angle_Err[3],Angle_Err[4],Angle_Err[5]);
    printf("Angular_Err:%4f,%4f,%4f,%4f,%4f,%4f\n",Angular_Err[0],Angular_Err[1],Angular_Err[2],Angular_Err[3],Angular_Err[4],Angular_Err[5]);
   Counter[0] = t/0.005;
   printf("Counter:%4f\n",Counter[0]);
}



#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*======================================================*
 * See sfuntmpl_doc.c for the optional S-function methods *
 *======================================================*/

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
