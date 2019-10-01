#include "GAMER.h"
#include "TestProb.h"   // for MODEL_MSTAR

// declare as static so that other functions cannot invoke them directly and must use the function pointers
static bool Flu_ResetByUser_Func( real fluid[], const double X, const double Y, const double Z, const double Time,
                                  const int lv, double AuxArray[] );
static void Flu_ResetByUser_API( const int lv, const int FluSg, const double TTime );

// these function pointers may be overwritten by various test problem initializers
bool (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double X, const double Y, const double Z, const double Time,
                                  const int lv, double AuxArray[] ) = Flu_ResetByUser_Func;
void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime ) = Flu_ResetByUser_API;


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
bool Flu_ResetByUser_Func( real fluid[], const double X, const double Y, const double Z, const double Time,
                           const int lv, double AuxArray[] )
{

// Example : reset fluid variables to extremely small values if the cell is within a specific sphere
   /*
   const double dr[3]     = { X - amr->BoxCenter[0],
                              Y - amr->BoxCenter[1],
                              Z - amr->BoxCenter[2] };
   const real r           = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   const real TRad        = 0.3;
   const real MaxDens     = 1.0e15;
   const real MaxPres     = 1.0e15;

   if ( r <= TRad )
   {
//    set active scalars
      fluid[DENS] = MaxDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = MaxPres / ( GAMMA-(real)1.0 );

//    set passive scalars

      return true;
   }

   else
      return false;
   */

   return false;

} // FUNCTION : Flu_ResetByUser_Func



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_API
// Description :  API for resetting the fluid array
//
// Note        :  1. Enabled by the runtime option "OPT__RESET_FLUID"
//                2. Invoked by either "Flu_AdvanceDt()" or "Gra_AdvanceDt()" using the function pointer
//                   "Flu_ResetByUser_API_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                3. Currently NOT applied to the input uniform array
//                   --> Init_ByFile() does NOT call this function
//                4. Currently does not work with "OPT__OVERLAP_MPI"
//                5. The function pointer "Flu_ResetByUser_Func_Ptr" points to "Flu_ResetByUser_Func()" by default
//                   but may be overwritten by various test problem initializers
//
// Parameter   :  lv    : Target refinement level
//                FluSg : Target fluid sandglass
//                TTime : Target physical time
//-------------------------------------------------------------------------------------------------------
void Flu_ResetByUser_API( const int lv, const int FluSg, const double TTime )
{
#  ifdef MODEL_MSTAR
   if (TTime > Time2Accrete) {
    
      // 1.0 MPI_AllReduce (or MPI_Reduce) to distribute d_mstar to d_mstar_sum
      MPI_Allreduce(&d_MStar,   &d_MStar_SUM,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&d_Star_J,  &d_Star_J_SUM,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(d_Star_Mom, d_Star_Mom_SUM, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
      // 2.0 reset M, J, Mom and GM  
      M_STAR += d_MStar_SUM  ;
      STAR_J += d_Star_J_SUM ; 
      for (int n=0; n<3; n++) Star_Mom[n] += d_Star_Mom_SUM[n] ;
      
      ExtAcc_AuxArray[3] = (real) NEWTON_G * M_STAR;

      const double _M_Star = 1.0/M_STAR ;
      
      // 3.0 update star's location
      /*
      const double dt = dTime_AllLv[0];
      double star_pos_cyl[3] = { ExtAcc_AuxArray[0], ExtAcc_AuxArray[1], ExtAcc_AuxArray[2] } ;
      double star_pos_crt[3] ;
      Aux_Coord_Adopted2CartesianCoord(star_pos_cyl, star_pos_crt);
      
      // move the star based on d_Star_Mom_SUM
      for (int n=0; n<3; n++) star_pos_crt[n] += Star_Mom[n]*_M_Star * dt ;
      
      //
      Aux_Coord_Cartesian2AdoptedCoord(star_pos_crt, star_pos_cyl);
      for (int n=0; n<3; n++) {
         Star_Pos[n]        = star_pos_cyl[n] ;
         ExtAcc_AuxArray[n] = Star_Pos[n] ;
      }
      */
         
   }
#  endif // #ifdef MODEL_MSTAR
   
   
   /*
// check
   if ( Flu_ResetByUser_Func_Ptr == NULL )
   {
      OPT__RESET_FLUID = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : OPT__RESET_FLUID has been disabled since Flu_ResetByUser_Func_Ptr == NULL !!\n" );

      return;
   }
   */

   const double *dh      = amr->dh[lv];
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   const real   Gamma_m1 = GAMMA - (real)1.0;
   const real  _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif

   bool   Reset;
   real   fluid[NCOMP_TOTAL];
   double X, Y, Z, X0, Y0, Z0;


#  pragma omp parallel for private( Reset, fluid, X, Y, Z, X0, Y0, Z0 ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      X0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh[0];
      Y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh[1];
      Z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh[2];

      for (int k=0; k<PS1; k++)  {  Z = Z0 + k*dh[2];
      for (int j=0; j<PS1; j++)  {  Y = Y0 + j*dh[1];
      for (int i=0; i<PS1; i++)  {  X = X0 + i*dh[0];

         for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[FluSg][lv][PID]->fluid[v][k][j][i];

//       reset this cell
         Reset = Flu_ResetByUser_Func_Ptr( fluid, X, Y, Z, TTime, lv, NULL );

//       operations necessary only when this cell has been reset
         if ( Reset )
         {
#           if ( MODEL == HYDRO  ||  MODEL == MHD )
//          check minimum density and pressure
            fluid[DENS] = FMAX( fluid[DENS], (real)MIN_DENS );
            fluid[ENGY] = CPU_CheckMinPresInEngy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                                  Gamma_m1, _Gamma_m1, MIN_PRES );

//          calculate the dual-energy variable (entropy or internal energy)
#           if   ( DUAL_ENERGY == DE_ENPY )
            fluid[ENPY] = CPU_Fluid2Entropy( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], Gamma_m1 );
#           elif ( DUAL_ENERGY == DE_EINT )
#           error : DE_EINT is NOT supported yet !!
#           endif

//          floor and normalize passive scalars
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = FMAX( fluid[v], TINY_NUMBER );

            if ( OPT__NORMALIZE_PASSIVE )
               CPU_NormalizePassive( fluid[DENS], fluid+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#           endif
#           endif // if ( MODEL == HYDRO  ||  MODEL == MHD )

//          store the reset values
            for (int v=0; v<NCOMP_TOTAL; v++)   amr->patch[FluSg][lv][PID]->fluid[v][k][j][i] = fluid[v];
         } // if ( Reset )

      }}} // i,j,k
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   
   
} // FUNCTION : Flu_ResetByUser_API
