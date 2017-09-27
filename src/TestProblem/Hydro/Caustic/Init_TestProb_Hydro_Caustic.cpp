#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Caustic_VelPeak;      // peak velocity
static double Caustic_Dens;         // background density
static double Caustic_Pres;         // background pressure
static int    Caustic_Dir;          // spatial direction: (0/1) --> (x/diagonal)
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

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 1\" (i.e., periodic BC) !!\n" );

#  if ( COORDINATE == CARTESIAN )
   if (  !Mis_CompareRealValue( amr->dh[0][0], amr->dh[0][1], NULL, false )  ||
         !Mis_CompareRealValue( amr->dh[0][0], amr->dh[0][2], NULL, false )    )
      Aux_Error( ERROR_INFO, "only work with cubic cells (dh[lv=0] = (%20.14e, %20.14e, %20.14e)) !!",
                 amr->dh[0][0], amr->dh[0][1], amr->dh[0][2] );
#  else
      Aux_Error( ERROR_INFO, "only work with the Cartesian coordinates !!\n" );
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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Caustic_VelPeak",   &Caustic_VelPeak,       -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Caustic_Dens",      &Caustic_Dens,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Caustic_Pres",      &Caustic_Pres,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Caustic_Dir",       &Caustic_Dir,            0,            0,                1                 );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default values


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = ( Caustic_Dir == 0 ) ? 3.0 : 2.0;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if (  Caustic_Dir == 1  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
         Aux_Error( ERROR_INFO, "simulation domain must be cubic for %s = %d !!\n", "Caustic_Dir", Caustic_Dir );


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",     TESTPROB_ID     );
      Aux_Message( stdout, "  peak velocity       = %13.7e\n", Caustic_VelPeak );
      Aux_Message( stdout, "  background density  = %13.7e\n", Caustic_Dens    );
      Aux_Message( stdout, "  background pressure = %13.7e\n", Caustic_Pres    );
      Aux_Message( stdout, "  spatial direction   = %s\n",     ( Caustic_Dir == 0 ) ? "X" :
                                                               ( Caustic_Dir == 1 ) ? "Diagonal" :
                                                                                      "Unknown" );
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

// x direction
   if ( Caustic_Dir == 0 )
   {
      const double WaveK = 2.0*M_PI/amr->BoxSize[0];

      fluid[MOMX] = Caustic_Dens*Caustic_VelPeak*sin( WaveK*x );
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
   }

// diagonal direction
   else if ( Caustic_Dir == 1 )
   {
      const double WaveK = 2.0*M_PI/amr->BoxSize[0]*sqrt(3.0);
      const double r     = 1.0/sqrt(3.0)*( x + y + z );

      fluid[MOMX] = Caustic_Dens*Caustic_VelPeak*sin( WaveK*r ) / sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
   }

   else
      Aux_Error( ERROR_INFO, "Caustic_Dir = %d is NOT supported [0/1] !!\n", Caustic_Dir );

   fluid[DENS] = Caustic_Dens;
   fluid[ENGY] = Caustic_Pres/(GAMMA-1.0) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Caustic
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Caustic()
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

} // FUNCTION : Init_TestProb_Hydro_Caustic
