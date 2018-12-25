#include "GAMER.h"
#include "TestProb.h"

static void Init_ExternalAcc() ;
static void BC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
                
                
// problem-specific global variables
// =======================================================================================
static double GM ;      // set code unit
static double R_0;      // set code unit
static double const_R;  // set code unit

static double slope_q;
static double slope_p;
static double T_0;      // T_0 ~ H^2
static double rho_0;
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
   Aux_Error( ERROR_INFO, "VSI only wroks with Cylindrical Coordinate !!\n" );
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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   // default to a Keplerian disk
   ReadPara->Add( "GM",                &GM,                    1.0,         Eps_double,        NoMax_double      );
   ReadPara->Add( "R_0",               &R_0,                   1.0,         Eps_double,        NoMax_double      );
   ReadPara->Add( "const_R",           &const_R,               1.0,         Eps_double,        NoMax_double      );
   ReadPara->Add( "slope_q",           &slope_q,               -2.0,        Eps_double,        NoMax_double      );
   ReadPara->Add( "slope_p",           &slope_p,               -1.5,        Eps_double,        NoMax_double      );
   ReadPara->Add( "T_0",               &T_0,                   4e-2,        Eps_double,        NoMax_double      );
   ReadPara->Add( "rho_0",             &rho_0,                 1.0,         Eps_double,        NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 100.0;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  GM                        = %13.7e\n", GM          );
      Aux_Message( stdout, "  R_0                       = %13.7e\n", R_0         );
      Aux_Message( stdout, "  const_R                   = %13.7e\n", const_R     );
      Aux_Message( stdout, "  Temperature slope index   = %13.7e\n", slope_q     );
      Aux_Message( stdout, "  Density slope index       = %13.7e\n", slope_p     );
      Aux_Message( stdout, "  T_0                       = %13.7e\n", T_0         );
      Aux_Message( stdout, "  rho_0                     = %13.7e\n", rho_0       );
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

   const real R_norm      = x / R_0; 
   const real rho_mid     = rho_0 * POW(R_norm, slope_p) ;
   const real temperature = T_0 * POW(R_norm, slope_q) ;
   const real cs_square   = const_R * temperature ; 
   const real _sph_r      = 1.0 / SQRT(x*x + z*z);
   const real _R_norm     = 1.0 / R_norm ;
   const real rho         = rho_mid * EXP( GM/cs_square * (_sph_r - _R_norm) );
   const real omega_kep   = SQRT(GM*_sph_r);
   const real H           = SQRT(cs_square)/omega_kep;
   
   // set up random number gen
   /*
   RandomNumber_t *RNG = NULL;
   RNG = new RandomNumber_t( 1 );
   double RanVel = RNG->GetValue( 0, -1.0*vel_amp, vel_amp ) ;
   */

   fluid[DENS] = rho;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = rho* x*omega_kep * SQRT(1.0 + (slope_q+slope_p)*SQR(H/x) + slope_q*(1.0-x*_sph_r) ) ;
   fluid[MOMZ] = 0.0;
   
   //fluid[MOMY] +=  RanVel ;

} // FUNCTION : SetGridIC


#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc() {
   
   ExtAcc_AuxArray[0] = GM ;

} // FUNCTION : Init_ExternalAcc
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition to the analytical solution
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{
   // this BC is only useful for checking the corectness of IC
   const real R_norm      = x / R_0; 
   const real rho_mid     = rho_0 * POW(R_norm, slope_p) ;
   const real temperature = T_0 * POW(R_norm, slope_q) ;
   const real cs_square   = const_R * temperature ; 
   const real _sph_r      = 1.0 / SQRT(x*x + z*z);
   const real _R_norm     = 1.0 / R_norm ;
   const real rho         = rho_mid * EXP( GM/cs_square * (_sph_r - _R_norm) );
   const real omega_kep   = SQRT(GM*_sph_r);
   const real H           = SQRT(cs_square)/omega_kep;

   fluid[DENS] = rho;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = rho* x*omega_kep * SQRT(1.0 + (slope_q+slope_p)*SQR(H/x) + slope_q*(1.0-x*_sph_r) ) ;
   fluid[MOMZ] = 0.0;

} // FUNCTION : BC




#endif // #if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_VSI
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_VSI()
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
#ifdef GRAVITY
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;       // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
#endif //GRAVITY
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_VSI
