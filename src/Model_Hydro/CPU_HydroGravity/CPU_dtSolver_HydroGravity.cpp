#include "GAMER.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )




//-----------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver_HydroGravity
// Description :  Estimate the evolution time-step (dt) required for the hydro gravity solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the free-fall time of the maximum gravitational acceleration
//
// Parameter   :  dt_Array        : Array to store the minimum dt in each target patch
//                Pot_Array       : Array storing the prepared potential data of each target patch
//                Corner_Array    : Array storing the physical corner coordinates of each patch
//                NPatchGroup     : Number of target patch groups
//                dh              : Grid size
//                Safety          : dt safety factor
//                P5_Gradient     : Use 5-points stencil to evaluate the potential gradient
//                GravityType     : Types of gravity --> self-gravity, external gravity, both
//                ExtAcc_AuxArray : Auxiliary array for adding external acceleration
//                ExtAcc_Time     : Physical time for adding the external acceleration
//
// Return      :  dt_Array
//-----------------------------------------------------------------------------------------
void CPU_dtSolver_HydroGravity( real dt_Array[],
                                const real Pot_Array[][ CUBE(GRA_NXT) ],
                                const double Corner_Array[][3],
                                const int NPatchGroup, const real dh[], const real Safety, const bool P5_Gradient,
                                const OptGravityType_t GravityType, const double ExtAcc_AuxArray[],
                                const double ExtAcc_Time )
{

// check
#  ifdef GAMER_DEBUG
   if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  &&  ExtAcc_Time < 0.0 )
      Aux_Error( ERROR_INFO, "Incorrect ExtAcc_Time (%14.7e) !!\n", ExtAcc_Time );
#  endif // #ifdef GAMER_DEBUG


   const int  NPatch       = NPatchGroup*8;
   const real Const_8      = (real)8.0;
   const int  did1[3]      = { 1, GRA_NXT, SQR(GRA_NXT) };
   const int  did2[3]      = { 2*did1[0], 2*did1[1], 2*did1[2] };
   
   real Gra_Const[3] ;     // ### this should be const real ...
   if (P5_Gradient)  for (int d=0; d<3; d++) Gra_Const[d] = -1.0/(12.0*dh[d]) ; 
   else              for (int d=0; d<3; d++) Gra_Const[d] = -1.0/(2.0 *dh[d])  ; 

   real   Acc[3], dx_Acc_Min;
   double x, y, z;
   int    id;
   
   real dx[3] = {dh[0], dh[1], dh[2]} ;
   real geo_factor = (real) 1.0 ;
   
#  if (COORDINATE == CYLINDRICAL)
   real radius ;
#  endif


// loop over all patches
#  pragma omp parallel for private( Acc, dx_Acc_Min, x, y, z, id ) schedule( runtime )
   for (int P=0; P<NPatch; P++)
   {
      dx_Acc_Min = HUGE_NUMBER;

      for (int k=GRA_GHOST_SIZE, kk=0; k<GRA_NXT-GRA_GHOST_SIZE; k++, kk++)
      for (int j=GRA_GHOST_SIZE, jj=0; j<GRA_NXT-GRA_GHOST_SIZE; j++, jj++)
      for (int i=GRA_GHOST_SIZE, ii=0; i<GRA_NXT-GRA_GHOST_SIZE; i++, ii++)
      {
         id     = ( k*GRA_NXT + j )*GRA_NXT + i;
         Acc[0] = (real)0.0;
         Acc[1] = (real)0.0;
         Acc[2] = (real)0.0;
         
#        if (COORDINATE == CYLINDRICAL)
         radius     = Corner_Array[P][0] + (double)(ii*dh[0]);
         geo_factor = (real)1.0 / radius; 
         dx[1]      = radius * dh[1] ;
#        endif

//       external gravity
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         {
            x = Corner_Array[P][0] + (double)ii*dh[0];
            y = Corner_Array[P][1] + (double)jj*dh[1];
            z = Corner_Array[P][2] + (double)kk*dh[2];

            CPU_ExternalAcc( Acc, x, y, z, ExtAcc_Time, ExtAcc_AuxArray );
         }

//       self-gravity
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
         {
            if ( P5_Gradient )
            {
               Acc[0] += Gra_Const[0] * ( -         Pot_Array[P][ id + did2[0] ] +         Pot_Array[P][ id - did2[0] ]
                                          + Const_8*Pot_Array[P][ id + did1[0] ] - Const_8*Pot_Array[P][ id - did1[0] ] );
               Acc[1] += Gra_Const[1] * ( -         Pot_Array[P][ id + did2[1] ] +         Pot_Array[P][ id - did2[1] ]
                                          + Const_8*Pot_Array[P][ id + did1[1] ] - Const_8*Pot_Array[P][ id - did1[1] ] ) * geo_factor;
               Acc[2] += Gra_Const[2] * ( -         Pot_Array[P][ id + did2[2] ] +         Pot_Array[P][ id - did2[2] ] 
                                          + Const_8*Pot_Array[P][ id + did1[2] ] - Const_8*Pot_Array[P][ id - did1[2] ] );
            }

            else
            {
               Acc[0] += Gra_Const[0] * ( Pot_Array[P][ id + did1[0] ] - Pot_Array[P][ id - did1[0] ] );
               Acc[1] += Gra_Const[1] * ( Pot_Array[P][ id + did1[1] ] - Pot_Array[P][ id - did1[1] ] ) * geo_factor ;
               Acc[2] += Gra_Const[2] * ( Pot_Array[P][ id + did1[2] ] - Pot_Array[P][ id - did1[2] ] );
            }
         } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )

//       get the maximum acceleration
         for (int d=0; d<3; d++)    dx_Acc_Min = FMIN( dx_Acc_Min, FABS(dx[d]/Acc[d]) );
      } // i,j,k

//    get the minimum dt
      dt_Array[P] = Safety*SQRT( (real)2.0*dx_Acc_Min );

   } // for (int P=0; P<NPatch; P++)

} // FUNCTION : CPU_dtSolver_HydroGravity



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
