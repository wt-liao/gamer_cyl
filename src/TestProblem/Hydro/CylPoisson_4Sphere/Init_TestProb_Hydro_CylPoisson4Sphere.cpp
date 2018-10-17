#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Sphere_Dens;          // sphere density
static double Sphere_Radius;
// =======================================================================================


static void Output_L1Error_CylPoisson4Sphere();

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
   ReadPara->Add( "Sphere_Dens",       &Sphere_Dens,           10.0,          Eps_double,       NoMax_double      );
//   ReadPara->Add( "Sphere_Radius",     &Sphere_Radius,         -1.0,          Eps_double,       NoMax_double      );
//   ReadPara->Add( "Sphere_Center_X",   &Sphere_Center[0],      NoDef_double,  NoMin_double,     NoMax_double      );
//   ReadPara->Add( "Sphere_Center_Y",   &Sphere_Center[1],      NoDef_double,  NoMin_double,     NoMax_double      );
//   ReadPara->Add( "Sphere_Center_Z",   &Sphere_Center[2],      NoDef_double,  NoMin_double,     NoMax_double      );


   ReadPara->Read( FileName );

   delete ReadPara;

// set the default 
   Sphere_Radius = 3.5;


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
      Aux_Message( stdout, "  sphere mass density       = %13.7e\n", Sphere_Dens  );
      Aux_Message( stdout, "  sphere radius             = %13.7e\n", Sphere_Radius);
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
   const double r_center = amr->BoxCenter[0] ;
   const double z_center = amr->BoxCenter[2] ;
   const double x1_cyl[4] = {r_center, r_center, r_center, r_center} ;
   const double x2_cyl[4] = {0, 0.5 * amr->BoxCenter[1], amr->BoxCenter[1], 1.5*amr->BoxCenter[1]} ;
   const double x3_cyl[4] = {z_center, z_center, z_center, z_center} ;
   
   
   const double x1 = x * COS(y) ;
   const double x2 = x * SIN(y) ;
   const double x3 = z ;
   
   double r, Sphere_C1, Sphere_C2, Sphere_C3 ;
   
   fluid[ENGY] = background_engy;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[DENS] = background_dens;
   
   for (int count=0; count<4; count++) {
      Sphere_C1 = x1_cyl[count] * COS(x2_cyl[count]);
      Sphere_C2 = x1_cyl[count] * SIN(x2_cyl[count]);
      Sphere_C3 = x3_cyl[count];
      
      r = SQRT( SQR(x1-Sphere_C1) + SQR(x2-Sphere_C2) + SQR(x3-Sphere_C3) );
      
      if (r <= Sphere_Radius) fluid[DENS] += Sphere_Dens;
   }

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    : Output_L1Error_CylPoisson
// Description : 
//-------------------------------------------------------------------------------------------------------
void Output_L1Error_CylPoisson4Sphere(){
#ifdef GRAVITY
   double analytic_sum = 0.0, L1_error = 0.0, total_error = 0.0;
   double global_analytic_sum = 0.0, global_total_error = 0.0;
   double L1_local = 0.0, L1_global = 0.0;
   double analytic_sol ;
   
   
   int FluSg, PotSg;
   double x_cyl, y_cyl, z_cyl, x_crt, y_crt, z_crt, dist2center, dV;
      
   /// set up for sphere coordinate
   const double r_center = amr->BoxCenter[0] ;
   const double z_center = amr->BoxCenter[2] ;
   const double x1_cyl[4] = {r_center, r_center, r_center, r_center} ;
   const double x2_cyl[4] = {0, 0.5 * amr->BoxCenter[1], amr->BoxCenter[1], 1.5*amr->BoxCenter[1]} ;
   const double x3_cyl[4] = {z_center, z_center, z_center, z_center} ;
   // sphere in crt coordinate 
   double sphere_x1[4], sphere_x2[4], sphere_x3[4];
   
   // sphere center in cartesian 
   for (int count=0; count<4; count++) {
      sphere_x1[count] = x1_cyl[count] * COS(x2_cyl[count]);
      sphere_x2[count] = x1_cyl[count] * SIN(x2_cyl[count]);
      sphere_x3[count] = x3_cyl[count];
   }
   
   
   // ## currently only for no amr
   for (int lv=0; lv<NLEVEL; lv++) {
      
      const double *dh = amr->dh[lv];
      FluSg = amr->FluSg[lv];
      PotSg = amr->PotSg[lv];
      
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++) {
         if ( amr->patch[0][lv][PID]->son == -1 ) {
               
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++) 
            {
               x_cyl = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 0, i ) ;
               y_cyl = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 1, j ) ;
               z_cyl = Aux_Coord_CellIdx2AdoptedCoord( lv, PID, 2, k ) ;
               dV    = Aux_Coord_CellIdx2Volume(lv, PID, i, j, k) ;
                  
               x_crt = x_cyl * COS(y_cyl) ;
               y_crt = x_cyl * SIN(y_cyl) ;
               z_crt = z_cyl ;
               
               analytic_sol = 0.0 ;
                              
               for (int count=0; count<4; count++) {
                  dist2center = SQRT( SQR(x_crt-sphere_x1[count]) + SQR(y_crt-sphere_x2[count]) + SQR(z_crt-sphere_x3[count]) );
                  
                  if ( dist2center <= Sphere_Radius ) 
                     analytic_sol += -2.0/3.0*M_PI*Sphere_Dens*( 3.0*SQR(Sphere_Radius) - SQR(dist2center) ) ;
                  else 
                     analytic_sol += -4.0/3.0*M_PI*Sphere_Dens*CUBE(Sphere_Radius)/dist2center ;
               }
               
                  
               analytic_sum += analytic_sol ;
               total_error += FABS( amr->patch[PotSg][lv][PID]->pot[k][j][i] - analytic_sol ) ;
               L1_local    += FABS( (amr->patch[PotSg][lv][PID]->pot[k][j][i] - analytic_sol)/analytic_sol*dV );
                  
            }// for k, j, i
         } // if leaf node
      } // for PID0; PID
      
   } // for lv
   
   //### MPI_Reduce: total_error
   //### MPI_reduce: analytic_sum
   
   MPI_Reduce(&total_error , &global_total_error , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&analytic_sum, &global_analytic_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&L1_local    , &L1_global          , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   
   if (MPI_Rank == 0)   {
      L1_error = global_total_error / FABS(global_analytic_sum) ;
      Aux_Message( stdout, "[Cylindrical Poisson Solver] L1 error (old version) = %20.12f. \n", L1_error );
      Aux_Message( stdout, "[Cylindrical Poisson Solver] L1 error (new version) = %20.12f. \n", L1_global );
   }
#endif // GRAVITY     
}


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
void Init_TestProb_Hydro_CylPoisson4Sphere()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
#ifdef GRAVITY
   Output_User_Ptr          = Output_L1Error_CylPoisson4Sphere;
#endif
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave
