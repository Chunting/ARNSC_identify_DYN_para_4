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
#define pi 3.14159265

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
 * ssSetNumSampleTimes(S, 1):  Specify the number of sample times that an S-Function block has.
 *
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
    
    if (!ssSetNumInputPorts(S, 16)) return;
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
    
    ssSetInputPortWidth(S, 6, 6);
    ssSetInputPortDirectFeedThrough(S, 6, 1);
    ssSetInputPortWidth(S, 7, 3);
    ssSetInputPortDirectFeedThrough(S, 7, 1);

    ssSetInputPortWidth(S, 8, 6);
    ssSetInputPortDirectFeedThrough(S, 8, 1);
    
    ssSetInputPortWidth(S, 9, 3);
    ssSetInputPortDirectFeedThrough(S, 9, 1);
    ssSetInputPortWidth(S, 10, 3);
    ssSetInputPortDirectFeedThrough(S, 10, 1);
    ssSetInputPortWidth(S, 11, 4);
    ssSetInputPortDirectFeedThrough(S, 11, 1);
    ssSetInputPortWidth(S, 12, 1);
    ssSetInputPortDirectFeedThrough(S, 12, 1);
    ssSetInputPortWidth(S, 13, 5);
    ssSetInputPortDirectFeedThrough(S, 13, 1);
    ssSetInputPortWidth(S, 14, 6);
    ssSetInputPortDirectFeedThrough(S, 14, 1);
    ssSetInputPortWidth(S, 15, 13);
    ssSetInputPortDirectFeedThrough(S, 15, 1);
    
    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
     */
    if (!ssSetNumOutputPorts(S, 6)) return; //设置输出变量的个数
    
    ssSetOutputPortWidth(S, 0, 1);
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);
    ssSetOutputPortWidth(S, 4, 1);
    ssSetOutputPortWidth(S, 5, 3);
    
    
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
    static real_T fai = pi/4.0;
    
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
    real_T theta0;
    real_T theta0_U;
    real_T theta0_D;
    real_T theta1_U;
    real_T theta2_U;
    real_T theta3_U;
    real_T theta1_U_t_1;
    real_T theta2_U_t_1;
    real_T theta3_U_t_1;
    real_T thetaE_U;
    real_T theta1_D;
    real_T theta2_D;
    real_T theta3_D;
    real_T theta1_D_t_1;
    real_T theta2_D_t_1;
    real_T theta3_D_t_1;
    
    real_T dtheta0;
    real_T dtheta0_U;
    real_T dtheta0_D;
    real_T dtheta1_U;
    real_T dtheta2_U;
    real_T dtheta3_U;
    real_T dtheta1_D;
    real_T dtheta2_D;
    real_T dtheta3_D;

    real_T L_m[1][1];

    time_T t = ssGetT(S);

    /*********** InputRealPtrsType ssGetInputPortRealSignalPtrs(SimStruct *S,int_T port) ***************/
    /***
     * S： SimStruct representing an S-Function block.
     * port: Index of the port whose signal is required.
     ***/
    // Get pointers to signals of type double connected to an input port.
    InputRealPtrsType 	 Base_x_in = ssGetInputPortRealSignalPtrs( S, 0 );  //基座X方向位移
    InputRealPtrsType 	 Base_y_in = ssGetInputPortRealSignalPtrs( S, 1 ); //基座Y方向位移
    InputRealPtrsType 	 Theta0_in = ssGetInputPortRealSignalPtrs( S, 2 ); //基座绕Z轴转动
    InputRealPtrsType 	 Base_vx_in = ssGetInputPortRealSignalPtrs( S, 3 ); //基座X方向速度
    InputRealPtrsType 	 Base_vy_in = ssGetInputPortRealSignalPtrs( S, 4 ); //基座Y方向速度
    InputRealPtrsType 	 dTheta0_in = ssGetInputPortRealSignalPtrs( S, 5 ); //基座绕Z轴转动角速度
    InputRealPtrsType 	 Theta_U_in = ssGetInputPortRealSignalPtrs( S, 6 ); // U杆三个关节的关节角度
    InputRealPtrsType 	 dTheta_U_in = ssGetInputPortRealSignalPtrs( S, 7 ); // U杆三个关节的关节角速度
    //  InputRealPtrsType 	 ddTheta_U_in = ssGetInputPortRealSignalPtrs( S, 8 ); // U杆三个关节的关节角加速度
    InputRealPtrsType 	 Theta_D_in = ssGetInputPortRealSignalPtrs( S, 8 );   // D杆三个关节的关节角度
    InputRealPtrsType 	 dTheta_D_in = ssGetInputPortRealSignalPtrs( S, 9 ); // D杆三个关节的关节角速度
    //   InputRealPtrsType 	 ddTheta_D_in = ssGetInputPortRealSignalPtrs( S, 11 ); // D杆三个关节的关节角加速度
    //   InputRealPtrsType 	 Link_CG_R_U_in = ssGetInputPortRealSignalPtrs( S, 12 );
    InputRealPtrsType 	 Pe_The_t_1_in = ssGetInputPortRealSignalPtrs( S, 10 );
    InputRealPtrsType 	 wR_U_t_1_in = ssGetInputPortRealSignalPtrs( S, 11 );  // pexy_ws_wc_out
    InputRealPtrsType 	 dTheta_T_in = ssGetInputPortRealSignalPtrs( S, 12 );  // Target 的角速度
    InputRealPtrsType 	 P_L_Mon_t_1 = ssGetInputPortRealSignalPtrs( S, 13 );  //Output (0-4), t-1时刻基座位姿和U臂的角度 角速度
    InputRealPtrsType 	 EndEffect_U_ve_w_in = ssGetInputPortRealSignalPtrs( S, 14 ); // 末端速度和角速度
    InputRealPtrsType 	 Iteration_in = ssGetInputPortRealSignalPtrs( S, 15 );
    
    
    // Get a pointer to an output signal of type double (real_T)
    real_T * Px_out = ssGetOutputPortRealSignal( S, 0 );
    real_T * Py_out = ssGetOutputPortRealSignal( S, 1 );
    real_T * L_out = ssGetOutputPortRealSignal( S, 2 );
    real_T * dTheta_Ue_out = ssGetOutputPortRealSignal( S, 3 );  // 末端角速度
    real_T * Pxs_Pyc_out = ssGetOutputPortRealSignal( S, 4 );
    real_T * Pe_The_out = ssGetOutputPortRealSignal( S, 5 );   // 末端位置和角度
  
    
    theta0 = * Theta0_in[0];
    
    theta0_U = theta0;
    theta0_D = theta0;
    
    theta1_U = * Theta_U_in[0];
    theta2_U = * Theta_U_in[1];
    theta3_U = * Theta_U_in[2];
    theta1_U_t_1 = * Theta_U_in[3];  //上个周期的实际角度
    theta2_U_t_1 = * Theta_U_in[4];
    theta3_U_t_1 = * Theta_U_in[5];
    
    
    theta1_D = * Theta_D_in[0];
    theta2_D = * Theta_D_in[1];
    theta3_D = * Theta_D_in[2];
    theta1_D_t_1 = * Theta_D_in[3];
    theta2_D_t_1 = * Theta_D_in[4];
    theta3_D_t_1 = * Theta_D_in[5];
    
    dtheta0 = * dTheta0_in[0];
    
    dtheta1_U = * dTheta_U_in[0];
    dtheta2_U = * dTheta_U_in[1];
    dtheta3_U = * dTheta_U_in[2];
    dtheta1_D = * dTheta_D_in[0];
    dtheta2_D = * dTheta_D_in[1];
    dtheta3_D = * dTheta_D_in[2];
    r0_0[0][0] = * Base_x_in[0];
    r0_0[1][0] = * Base_y_in[0];
    v0_0[0][0] = * Base_vx_in[0];
    v0_0[1][0] = * Base_vy_in[0];
    dTheta_U_RL_t_1[0] = * Iteration_in[0];
    dTheta_U_RL_t_1[1] = * Iteration_in[1];
    dTheta_U_RL_t_1[2] = * Iteration_in[2];
    dTheta_U_RL_t_1[3] = * Iteration_in[3];
    dTheta_U_RL_t_1[4] = * Iteration_in[4];
    dTheta_U_RL_t_1[5] = * Iteration_in[5];
    Q_t_1[0][0] = * Iteration_in[6];
    W_t_1[0][0] = * Iteration_in[7];
    W_t_1[0][1] = * Iteration_in[8];
    W_t_1[0][2] = * Iteration_in[9];
    W_t_1[0][3] = * Iteration_in[10];
    W_t_1[0][4] = * Iteration_in[11];
    W_t_1[0][5] = * Iteration_in[12];
    dtheta_T = * dTheta_T_in[0];
    ve_U[0][0]= * EndEffect_U_ve_w_in[0];  // 末端 X 方向速度
    ve_U[1][0]= * EndEffect_U_ve_w_in[1];  // 末端 Y 方向速度
    ve_U[2][0]= * EndEffect_U_ve_w_in[5];  // 末端绕 Z 方向角速度
    Px_t_1 =  * P_L_Mon_t_1[0];  // t-1时刻基座X方向位置
    Py_t_1 =  * P_L_Mon_t_1[1];  // t-1时刻基座Y方向位置
    L_t_1 =  * P_L_Mon_t_1[2];   // t-1时刻基座绕Z方向旋转角度
    dthetaE_t_1 =  * P_L_Mon_t_1[3];
    PcR_t_1 =  * P_L_Mon_t_1[4];
    ve_U_t_1[0][0] =  * wR_U_t_1_in[0];
    ve_U_t_1[1][0] =  * wR_U_t_1_in[1];
    dPIA_01 =  * wR_U_t_1_in[2];
    dPIA_11 =  * wR_U_t_1_in[3];
    re_the_U_t_1[0][0] =  * Pe_The_t_1_in[0];
    re_the_U_t_1[1][0] =  * Pe_The_t_1_in[1];
    re_the_U_t_1[2][0] =  * Pe_The_t_1_in[2];
    
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
    

//将线动量方程与角动量方程合并后，可表示为H_U*dTheta_U + H_D*dTheta_L + H_th0*dTheta0 = 0
//H_U、H_D为1X3矩阵
    H_U[0][0] = H_LLM_U[0][0] + H_PUM[0][0]*r_g[1][0] - H_PUM[1][0]*r_g[0][0];
    H_U[0][1] = H_LLM_U[0][1] + H_PUM[0][1]*r_g[1][0] - H_PUM[1][1]*r_g[0][0];
    H_U[0][2] = H_LLM_U[0][2] + H_PUM[0][2]*r_g[1][0] - H_PUM[1][2]*r_g[0][0];
    
    H_D[0][0] = H_LLM_D[0][0] + H_PDM[0][0]*r_g[1][0] - H_PDM[1][0]*r_g[0][0];
    H_D[0][1] = H_LLM_D[0][1] + H_PDM[0][1]*r_g[1][0] - H_PDM[1][1]*r_g[0][0];
    H_D[0][2] = H_LLM_D[0][2] + H_PDM[0][2]*r_g[1][0] - H_PDM[1][2]*r_g[0][0];
    
    H_th0 = H_LBth0 - (H_PBthy*r_g[0][0]-H_PBthx*r_g[1][0]);
    
    // 系统总体角动量？？？？？
    L_am = H_th0*dtheta0 + H_U[0][0]*dtheta1_U + H_U[0][1]*dtheta2_U + H_U[0][2]*dtheta3_U + H_D[0][0]*dtheta1_D + H_D[0][1]*dtheta2_D+H_D[0][2]*dtheta3_D;
    //L = H_LB[0][0]*v0_0[0][0]+H_LB[0][1]*v0_0[1][0]+ H_LBth0*dtheta0+H_LLM_U[0][0]*dtheta1_U+H_LLM_U[0][1]*dtheta2_U+H_LLM_U[0][2]*dtheta3_U+H_LLM_D[0][0]*dtheta1_D+H_LLM_D[0][1]*dtheta2_D+H_LLM_D[0][2]*dtheta3_D + (Px)*re_U[1][0]-(Py)*re_U[0][0];
    // 角动量
    L = H_LB[0][0]*v0_0[0][0] + H_LB[0][1]*v0_0[1][0] + H_LBth0*dtheta0 + H_LLM_U[0][0]*dtheta1_U + H_LLM_U[0][1]*dtheta2_U + H_LLM_U[0][2]*dtheta3_U + H_LLM_D[0][0]*dtheta1_D + H_LLM_D[0][1]*dtheta2_D + H_LLM_D[0][2]*dtheta3_D;
    printf("L:%4f\n",H_LB[0][0]*v0_0[0][0]+H_LB[0][1]*v0_0[1][0]+ H_LBth0*dtheta0+H_LLM_U[0][0]*dtheta1_U+H_LLM_U[0][1]*dtheta2_U+H_LLM_U[0][2]*dtheta3_U+H_LLM_D[0][0]*dtheta1_D+H_LLM_D[0][1]*dtheta2_D+H_LLM_D[0][2]*dtheta3_D);

    Px_out[0] =  Px;
    Py_out[0] =  Py;
    L_out[0] = L; // 系统角动量（不含捕获目标）
    dTheta_Ue_out[0] = ve_U[2][0] ;//dtheta0 + dtheta1_U + dtheta2_U + dtheta3_U;
    Pxs_Pyc_out[0] = (Px)*(sin(thetaE_U))-(Py)*(cos(thetaE_U));
    pexy_ws_wc_out[0] =  ve_U[0][0];  // 末端X位置
    pexy_ws_wc_out[1] =  ve_U[1][0];  // 末端Y位置
    pexy_ws_wc_out[2] =  -1.0*sin(thetaE_U)*ve_U[2][0];
    pexy_ws_wc_out[3] =  cos(thetaE_U)*ve_U[2][0];
    
    Pe_The_out[0] = re_U[0][0];
    Pe_The_out[1] = re_U[1][0];
    Pe_The_out[2] = thetaE_U;
    
    printf("omega_UD_dst:%4f,%4f,%4f,%4f,%4f,%4f\n",omega_UD_dst[0],omega_UD_dst[1],omega_UD_dst[2],omega_UD_dst[3],omega_UD_dst[4],omega_UD_dst[5]);
    printf("\n***********************************************\n");

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
