#ifdef GRAVITY

#include "GAMER.h"
#include "TestProb.h"


static void Init_ExternalAcc() ;
static bool Flu_ResetByUser( real fluid[], const double X, const double Y, const double Z, const double Time,
                             const int lv, double AuxArray[] ) ;
static void Aux_Record_User() ;

// problem-specific global variables
// =======================================================================================
#ifndef MODEL_MSTAR
static double M_STAR;
#endif
// =======================================================================================



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  if ( COORDINATE != CYLINDRICAL )
   Aux_Error( ERROR_INFO, "This test problem is currently only in Cylindrical coordinate !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__GRAVITY_TYPE != GRAVITY_SELF )
      Aux_Message( stderr, "WARNING : OPT__GRAVITY_TYPE != GRAVITY_SELF ??\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format (some handy constants are defined in TestProb.h):
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
#  ifndef MODEL_MSTAR
   ReadPara->Add( "M_STAR",            &M_STAR,                2.0e9,         Eps_double,       NoMax_double      );
#  endif
//   ReadPara->Add( "Sphere_Radius",     &Sphere_Radius,         -1.0,          Eps_double,       NoMax_double      );


   ReadPara->Read( FileName );

   delete ReadPara;

// set the default


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 7.0e-2;
   const long   End_Step_Default = 0 ;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

/*
   if ( !OPT__INIT_RESTRICT ) {
      OPT__INIT_RESTRICT = true;
      PRINT_WARNING( "OPT__INIT_RESTRICT", OPT__INIT_RESTRICT, FORMAT_BOOL );
   }
*/


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID  );
      Aux_Message( stdout, "  star mass                 = %13.7e\n", M_STAR       );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
   Aux_Message(stdout, "Initial condition read-in from UM_IC. \n");

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Func
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates in the adopted coordinate system
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser( real fluid[], const double X, const double Y, const double Z, const double Time,
                      const int lv, double AuxArray[] )
{
   // for MODEL_MSTAR, reset updates the M_star. Thus, no need to loop through XYZ;  
   // Put the reset directly in Flu_ResetByUser_API; 
   
#  ifdef SET_T_LIMIT_POPIII
   const double time_unit       = 5.02280842159966e+09;
   const double length_unit     = 1.49597870750767e+13;
   const double T_CMB           = 50 ;
   const double T_upper         = 5e4;
   const double m_ave_cgs       = Const_mH * (0.76 + 0.24*4) ;
   const double R               = (Const_kB/m_ave_cgs) * SQR(time_unit/length_unit) ;
   const double Gamma_m1        = GAMMA - 1.0;
   const double _Gamma_m1       = 1.0/Gamma_m1;
   const bool   CheckMinPres_No = false;
   
   double pres, pres_reset, ie_reset, T, dens, _dens; 
   
   dens  = fluid[DENS]; 
   _dens = 1 / dens ;
   pres  = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], 
                            Gamma_m1, CheckMinPres_No, MIN_PRES ); 
   T = pres*_dens/R;
   
   // reset for T_CMB
   if (T < T_CMB) {
      pres_reset  = dens*R*T_CMB;
      ie_reset    = pres_reset * _Gamma_m1;
      fluid[ENGY] = 0.5*(SQR(fluid[MOMX])+SQR(fluid[MOMY])+SQR(fluid[MOMZ]))*_dens + ie_reset;
      return true;
   }
   
   // reset for T_upper
   else if (T > T_upper) {
      pres_reset  = dens*R*T_upper;
      ie_reset    = pres_reset * _Gamma_m1;
      fluid[ENGY] = 0.5*(SQR(fluid[MOMX])+SQR(fluid[MOMY])+SQR(fluid[MOMZ]))*_dens + ie_reset;
      return true;
   }
   
#  endif // SET_T_LIMIT
   
#  ifdef SET_T_LIMIT_SG
   const double time_unit       = 5.02280842159966e+06;
   const double length_unit     = 1.49597870750767e+13;
   const double T_lower         = 10 ;
   const double T_upper         = 1e6;
   const double m_ave_cgs       = Const_mH ;
   const double R               = (Const_kB/m_ave_cgs) * SQR(time_unit/length_unit) ;
   const double Gamma_m1        = GAMMA - 1.0;
   const double _Gamma_m1       = 1.0/Gamma_m1;
   const bool   CheckMinPres_No = false;
   
   double pres, pres_reset, ie_reset, T, dens, _dens; 
   
   dens  = fluid[DENS]; 
   _dens = 1 / dens ;
   pres  = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], 
                            Gamma_m1, CheckMinPres_No, MIN_PRES ); 
   T = pres*_dens/R;
   
   // reset for T_lower
   if (T < T_lower) {
      pres_reset  = dens*R*T_lower;
      ie_reset    = pres_reset * _Gamma_m1;
      fluid[ENGY] = 0.5*(SQR(fluid[MOMX])+SQR(fluid[MOMY])+SQR(fluid[MOMZ]))*_dens + ie_reset;
      return true;
   }
   
   // reset for T_upper
   else if (T > T_upper) {
      pres_reset  = dens*R*T_upper;
      ie_reset    = pres_reset * _Gamma_m1;
      fluid[ENGY] = 0.5*(SQR(fluid[MOMX])+SQR(fluid[MOMY])+SQR(fluid[MOMZ]))*_dens + ie_reset;
      return true;
   }
   
#  endif // SET_T_LIMIT
   
   return false;
} // FUNCTION : Flu_ResetByUser_Func



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
// Description :  Initialize the external potential routines "CUPOT_ExternalAcc.cu / CPU_ExternalAcc.cpp"
//
// Note        :  Fill in the array "ExtAcc_AuxArray" here
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc()
{

// ExtAcc_AuxArray has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)

   ExtAcc_AuxArray[0] = (real) Star_Pos[0];
   ExtAcc_AuxArray[1] = (real) Star_Pos[1];
   ExtAcc_AuxArray[2] = (real) Star_Pos[2];
   ExtAcc_AuxArray[3] = (real) NEWTON_G * M_STAR;
   ExtAcc_AuxArray[4] = (real) 0.0;

} // FUNCTION : Init_ExternalAcc



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_User
// Description :  Record any user-specified information
//
// Note        :  1. Invoked by "main" using the function pointer "Aux_Record_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_User()
{
#  ifdef MODEL_MSTAR

   const char FileName[] = "Record__User";
   static bool FirstTime = true;

   if ( FirstTime )
   {
//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s%14s\n",  
                 "Time", "Step", "", "dt", "d_MStar", "GM", "Star_R", "Star_Theta", "Star_Z", "Star_J" );
         fclose( File_User );
      }

      FirstTime = false;
   }

// user-specified info
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e%14ld%3s%14.7e%14.7e%14.7e%14.6e%14.6e%14.6e%14.6e\n", 
               Time[0], Step, "", dTime_Base, d_MStar_SUM, ExtAcc_AuxArray[3], 
               ExtAcc_AuxArray[0], ExtAcc_AuxArray[1], ExtAcc_AuxArray[2], STAR_J );
               
      fclose( File_User );
   }
   
#  endif // MODEL_MSTAR
   
} // FUNCTION : Aux_Record_User


#  ifdef SUPPORT_GRACKLE
//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_GrackleOp
// Description :  Add the problem-specific fields
// 
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_GrackleOp()
{
#  if (defined GRACKLE_H2_SOBOLEV)
   Idx_alpha  = AddField( "Alpha",     NORMALIZE_NO );
   Idx_OpTauX = AddField( "Opacity_X", NORMALIZE_NO );
   Idx_OpTauY = AddField( "Opacity_Y", NORMALIZE_NO );
   Idx_OpTauZ = AddField( "Opacity_Z", NORMALIZE_NO );
   
#  elif (defined GRACKLE_H2_DISK)
   Idx_DiskTau = AddField( "DiskOpacity", NORMALIZE_NO );
   
#  endif
} // FUNCTION : AddNewField_GrackleOp

#endif // #ifdef SUPPORT_GRACKLE

#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_SGThinDisk()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = Aux_Record_User;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser;
   End_User_Ptr             = NULL;
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;       // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
   
#  if (defined SUPPORT_GRACKLE) && ((defined GRACKLE_H2_SOBOLEV) || (defined GRACKLE_H2_DISK))
   Init_Field_User_Ptr      = AddNewField_GrackleOp;
#  endif
   
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave

#endif // #ifdef GRAVITY
