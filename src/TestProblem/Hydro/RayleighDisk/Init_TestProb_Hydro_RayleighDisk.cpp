#include "GAMER.h"
#include "TestProb.h"

static void Init_ExternalAcc() ;
static void BC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
                
                
// problem-specific global variables
// =======================================================================================
static double Rayleigh_Slope;              // slope index = n; a = r^(-n)
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
   Aux_Error( ERROR_INFO, "Rayleigh Disk only wroks with Cylindrical Coordinate !!\n" );
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
   ReadPara->Add( "Rayleigh_Slope",    &Rayleigh_Slope,        2.0,         Eps_double,        NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 5.0e2;
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
      Aux_Message( stdout, "  Rayleigh disk slope index = %13.7e\n", Rayleigh_Slope );
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
   // Keplerian: n = 2
   // unstable disk: n > 3
   
   // g       = r^(-n)
   // v_theta = r^( (1-n)/2) )
   
   const real density  = (real) 50.0 ;
   const real n        = Rayleigh_Slope ;
   const real vel_amp  = (real) 0.0001 ;
   const real r_inner  = amr->BoxEdgeL[0] + 0.5*amr->dh[0][0];
   const real r_size   = amr->BoxEdgeR[0] - amr->BoxEdgeL[0] ;
   const real wave_num = 1.0 ;
   const real wave_k   = (2.0*M_PI) / r_size * wave_num ; 
   
   // set up random number gen
   RandomNumber_t *RNG = NULL;
   RNG = new RandomNumber_t( 1 );
   double RanVel = RNG->GetValue( 0, -1.0*vel_amp, vel_amp ) ;

   fluid[DENS] = density;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = density*POW( x, 0.5*(1-n) );
   fluid[MOMZ] = 0.0;
   
   //fluid[MOMY] +=  RanVel ;
   fluid[MOMY] += vel_amp*SIN( wave_k * (x-r_inner) ) ; 
   fluid[ENGY] = 0.5 * SQR(fluid[MOMY]) / density + density*0.0001 ;
   

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc() {
   
   ExtAcc_AuxArray[0] = Rayleigh_Slope ;

} // FUNCTION : Init_ExternalAcc




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

   const real density  = (real) 50.0 ;
   const real n        = Rayleigh_Slope ;
   
   fluid[DENS] = density;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = density*POW( x, 0.5*(1-n) );
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = 0.5 * SQR(fluid[MOMY]) / density + density*0.0001 ;

} // FUNCTION : BC




#endif // #if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_RayleighDisk
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_RayleighDisk()
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
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;       // option: OPT__GRAVITY_TYPE=2/3; example: SelfGravity/Init_ExternalAcc.cpp
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_RayleighDisk
