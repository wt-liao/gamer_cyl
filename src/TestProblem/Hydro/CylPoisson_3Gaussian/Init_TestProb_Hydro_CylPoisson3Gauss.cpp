#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
//static double Sphere_Radius;
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
//   ReadPara->Add( "Sphere_Dens",       &Sphere_Dens,           10.0,          Eps_double,       NoMax_double      );
//   ReadPara->Add( "Sphere_Radius",     &Sphere_Radius,         -1.0,          Eps_double,       NoMax_double      );
//   ReadPara->Add( "Sphere_Center_X",   &Sphere_Center[0],      NoDef_double,  NoMin_double,     NoMax_double      );
//   ReadPara->Add( "Sphere_Center_Y",   &Sphere_Center[1],      NoDef_double,  NoMin_double,     NoMax_double      );
//   ReadPara->Add( "Sphere_Center_Z",   &Sphere_Center[2],      NoDef_double,  NoMin_double,     NoMax_double      );


   ReadPara->Read( FileName );

   delete ReadPara;

// set the default 
   //Sphere_Radius = 3.5;


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
//      Aux_Message( stdout, "  sphere mass density       = %13.7e\n", Sphere_Dens  );
//      Aux_Message( stdout, "  sphere radius             = %13.7e\n", Sphere_Radius);
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
   
   const double background_dens = 0.0;
   const double background_engy = 0.0;
   
   /// set up for sphere coordinate
   const double r_center      = amr->BoxCenter[0] ;
   const double z_center      = amr->BoxCenter[2] ;
   const double x1_cyl[3]     = {1.0, 0.9, 1.0} ;
   const double x2_cyl[3]     = {0.0, 3.0*M_PI/4.0, 3.0*M_PI/2.0} ;
   const double x3_cyl[3]     = {z_center, z_center, z_center} ;
   const double sphere_den[3] = {2.0,  0.5,  1.0} ;
   const double sphere_sig[3] = {0.05, 0.05, 0.05} ;   
   
   const double x1 = x * COS(y) ;
   const double x2 = x * SIN(y) ;
   const double x3 = z ;
   
   double r, Sphere_C1, Sphere_C2, Sphere_C3 ;
   
   fluid[ENGY] = background_engy;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[DENS] = background_dens;
   
   for (int n=0; n<3; n++) {
      Sphere_C1 = x1_cyl[n] * COS(x2_cyl[n]);
      Sphere_C2 = x1_cyl[n] * SIN(x2_cyl[n]);
      Sphere_C3 = x3_cyl[n];
      
      r = SQRT( SQR(x1-Sphere_C1) + SQR(x2-Sphere_C2) + SQR(x3-Sphere_C3) );
      
      fluid[DENS] += sphere_den[n] * EXP( -SQR(r)/(2.0*SQR(sphere_sig[n])) ) / POW(2.0*M_PI*SQR(sphere_sig[n]), 1.5) ;
   }

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    : Output_L1Error_CylPoisson
// Description : 
//-------------------------------------------------------------------------------------------------------

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
void Init_TestProb_Hydro_CylPoisson3Gauss()
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
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave
