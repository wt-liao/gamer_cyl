#include "GAMER.h"
#include "CUFLU.h"

#ifdef MODEL_MSTAR
#include "TestProb.h"
#endif

#if (  !defined GPU  &&  MODEL == HYDRO  &&  \
       ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )

#if ( COORDINATE == CYLINDRICAL )
extern void GetCoord( const double Corner[], const real dh[], const int loop_size, real x_pos[], real face_pos[][2], 
                      const int i, const int j, const int k );
extern void GeometrySourceTerm( const real PriVar[], const real x_pos[], real GeoSource[] );
static void CurviFluxGrad( real dF[][NCOMP_TOTAL], const real x_pos[] );
static void GetFullStepGeoSource( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real* GeoSource, 
                                  const real dF[][NCOMP_TOTAL], const real* x_pos, const real* dt_dh2, const real dt_2, 
                                  const real Gamma_m1, const real MinPres, const int ID3 ) ;
extern void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                         const bool NormPassive, const int NNorm, const int NormIdx[],
                         const bool JeansMinPres, const real JeansMinPres_Coeff );
// for cooling
#ifdef COOLING
extern void CoolingFunc(real* cool_rate, const real PriVar[], const real x_pos[]);
#endif // COOLING
#endif

//### add MODEL_MSTAR in FullStepUpdate; 
//### -> this only works for CPU solver, since both d_MStar and Edge_x1_L are not passed into the solver

//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Parameter   :  Input            : Array storing the input initial data
//                Output           : Array to store the ouptut updated data
//                DE_Status        : Array to store the dual-energy status
//                Flux             : Array storing the input face-centered flux
//                                   --> Size is assumed to be N_FL_FLUX^3
//                dt               : Time interval to advance solution
//                dh               : Grid size
//                Gamma            : Ratio of specific heats
//                MinDens          : Minimum allowed density
//                MinPres          : Minimum allowed pressure
//                DualEnergySwitch : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive      : true --> normalize passive scalars so that the sum of their mass density
//                                            is equal to the gas mass density
//                NNorm            : Number of passive scalars to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx          : Target variable indices to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                         const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh[], const double Corner[],
                         const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                         const bool NormPassive, const int NNorm, const int NormIdx[] )
{

#  ifdef DUAL_ENERGY
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif
   const int  dID1[3]   = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh[3]  = { dt/dh[0], dt/dh[1], dt/dh[2] };

   int  ID1, ID2, ID3;
   real dF[3][NCOMP_TOTAL];
   
#  if ( COORDINATE == CYLINDRICAL )
   const real dt_2      = (real)0.5 * dt ;
   const real dt_dh2[3] = {dt_2/dh[0], dt_2/dh[1], dt_2/dh[2]};
   real x_pos[3], face_pos[1][2], GeoSource[NCOMP_TOTAL] ;
#  ifndef DUAL_ENERGY
   const real  Gamma_m1 = Gamma - (real)1.0; //### this has already been declared in DUAL_ENERGY
#  endif // #ifndef DUAL_ENERGY
#  ifdef MODEL_MSTAR
   const real Edge_x1_L = amr->BoxEdgeL[0];
   real dist2center, r_i;
   double d_star_mom_r, d_star_mom_theta, cos_theta, sin_theta;
   const double star_pos[3] = {ExtAcc_AuxArray[0], ExtAcc_AuxArray[1], ExtAcc_AuxArray[2]};
#  endif // MODEL_MSTAR
#  endif // COORDINATE == CYLINDRICAL

#  if ( NCOMP_PASSIVE > 0 )
   real Passive[NCOMP_PASSIVE];
#  endif


   for (int k1=0, k2=FLU_GHOST_SIZE;  k1<PS2;  k1++, k2++)
   for (int j1=0, j2=FLU_GHOST_SIZE;  j1<PS2;  j1++, j2++)
   for (int i1=0, i2=FLU_GHOST_SIZE;  i1<PS2;  i1++, i2++)
   {

      ID1 = (k1*N_FL_FLUX + j1)*N_FL_FLUX + i1;
      ID2 = (k1*PS2       + j1)*PS2       + i1;
      ID3 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)   dF[d][v] = Flux[ ID1+dID1[d] ][d][v] - Flux[ID1][d][v];
      
#     if (COORDINATE == CYLINDRICAL)
      GetCoord( Corner, dh, PS2, x_pos, face_pos, i1, j1, k1);
      CurviFluxGrad(dF, x_pos) ;
      GetFullStepGeoSource( Input, GeoSource, dF, x_pos, dt_dh2, dt_2, Gamma_m1, MinPres, ID3);
      
#     ifdef MODEL_MSTAR
      // only account for the flux from the inner most r-grid; be carful about ghost zone
      r_i         = x_pos[0]-0.5*dh[0] ;
      //dist2center = SQRT( SQR(r_i) + SQR(x_pos[2]) ) ;
      dist2center = SQRT( SQR(r_i) + SQR(star_pos[0]) - 2*r_i*star_pos[0]*cos(x_pos[1]-star_pos[1]) 
                        + SQR(x_pos[2]-star_pos[2]) ) ;
      
      if (x_pos[0] > Edge_x1_L && x_pos[0] < Edge_x1_L+dh[0] && dist2center < ACCRETE_RADIUS ) {
         //### Note that Flux = physical_flux*r_i
         d_MStar         += FMAX( -Flux[ID1][0][DENS], 0 ) * dt * (dh[1]*dh[2]) ; 
         
         if (Flux[ID1][0][DENS] < 0) {
            d_star_mom_r     = Flux[ID1][0][MOMX] * dt * (dh[1]*dh[2]) ;
            d_star_mom_theta = Flux[ID1][0][MOMY] * dt * (dh[1]*dh[2]) / r_i ;
            cos_theta        = cos(x_pos[1]);
            sin_theta        = sin(x_pos[1]);
         
            // momentum change in cartesian 
            d_Star_Mom[0]   += d_star_mom_r*cos_theta - d_star_mom_theta*sin_theta ;
            d_Star_Mom[1]   += d_star_mom_r*sin_theta + d_star_mom_theta*cos_theta ;
            d_Star_Mom[2]   += Flux[ID1][0][MOMZ] * dt * (dh[1]*dh[2]) ;
            d_Star_J        += d_star_mom_theta * r_i ;
         }
         
      } 
#     endif // MODEL_MSTAR
      
#     endif // COORDINATE == CYLINDRICAL

      for (int v=0; v<NCOMP_TOTAL; v++) {
         Output[v][ID2] = Input[v][ID3] - ( dF[0][v]*dt_dh[0] + dF[1][v]*dt_dh[1] + dF[2][v]*dt_dh[2] );
#        if (COORDINATE == CYLINDRICAL)
         Output[v][ID2] += GeoSource[v] * dt ;
#        endif
      }
      
      //### modified by wtl for debugging; only work for cyl coord
      /*
      double pres_debug1 = CPU_GetPressure( Input[DENS][ID3], Input[MOMX][ID3], Input[MOMY][ID3], 
                                            Input[MOMZ][ID3], Input[ENGY][ID3],
                                            Gamma_m1, false, NULL_REAL );
      
      double pres_debug  = CPU_GetPressure( Output[DENS][ID2], Output[MOMX][ID2], Output[MOMY][ID2], 
                                            Output[MOMZ][ID2], Output[ENGY][ID2],
                                            Gamma_m1, false, NULL_REAL );
                            
      if ( pres_debug < MIN_PRES ) 
         Aux_Message(stdout, "Weird prssure = %8.4e (before stepping: %8.4e) at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     pres_debug, pres_debug1, x_pos[0], x_pos[1], x_pos[2]);
                          
      if ( Output[DENS][ID2]<0 || !Aux_IsFinite(Output[DENS][ID2]) )
         Aux_Message(stdout, "Weird density = %8.4e (before stepping: %8.4e) at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     Output[DENS][ID2], Input[DENS][ID3], x_pos[0], x_pos[1], x_pos[2]);
                     
      if ( Output[ENGY][ID2]<0 || !Aux_IsFinite(Output[ENGY][ID2]) )
         Aux_Message(stdout, "Weird energy = %8.4e at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     Output[ENGY][ID2], x_pos[0], x_pos[1], x_pos[2]);
      
      if ( !Aux_IsFinite(Output[MOMX][ID2]) )
         Aux_Message(stdout, "Weird Mom_x = %8.4e at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     Output[MOMX][ID2], x_pos[0], x_pos[1], x_pos[2]);
                     
      if ( !Aux_IsFinite(Output[MOMY][ID2]) )
         Aux_Message(stdout, "Weird Mom_y = %8.4e at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     Output[MOMY][ID2], x_pos[0], x_pos[1], x_pos[2]);
                     
      if ( !Aux_IsFinite(Output[MOMZ][ID2]) )
         Aux_Message(stdout, "Weird Mom_z = %8.4e at (X,Y,Z)=(%6.3f, %6.3f, %6.3f). \n", 
                     Output[MOMZ][ID2], x_pos[0], x_pos[1], x_pos[2]);
      */
            

//    we no longer ensure positive density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//        because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
//    --> this consideration holds even when DUAL_ENERGY is adopted (e.g., when density is negative, even when DUAL_ENERGY is on,
//        we still want to try the 1st-order-flux correction before setting a floor value)
      /*
      Output[DENS][ID2] = FMAX( Output[DENS][ID2], MinDens );
      Output[ENGY][ID2] = CPU_CheckMinPresInEngy( Output[DENS][ID2], Output[MOMX][ID2], Output[MOMY][ID2], Output[MOMZ][ID2],
                                                  Output[ENGY][ID2], Gamma_m1, _Gamma_m1, MinPres );
      */


//    floor and normalize passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Output[v][ID2] = FMAX( Output[v][ID2], TINY_NUMBER );

      if ( NormPassive )
      {
         for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Output[ NCOMP_FLUID + v ][ID2];

         CPU_NormalizePassive( Output[DENS][ID2], Passive, NNorm, NormIdx );

         for (int v=0; v<NCOMP_PASSIVE; v++)    Output[ NCOMP_FLUID + v ][ID2] = Passive[v];
      }
#     endif


//    apply the dual-energy formalism to correct the internal energy
//    --> currently, even when UNSPLIT_GRAVITY is on (which would update the internal energy), we still invoke
//        CPU_DualEnergyFix() here and will fix the internal energy in the gravity solver for cells updated
//        by the dual-energy formalism (i.e., for cells with their dual-energy status marked as DE_UPDATED_BY_DUAL)
//    --> this feature might be modified in the future
#     ifdef DUAL_ENERGY
//    we no longer apply the minimum density and pressure checks here since we want to enable 1st-order-flux correction for that
      const bool CheckMinPres_No = false;
//    Output[DENS][ID2] = FMAX( Output[DENS][ID2], MinDens );

      CPU_DualEnergyFix( Output[DENS][ID2], Output[MOMX][ID2], Output[MOMY][ID2], Output[MOMZ][ID2],
                         Output[ENGY][ID2], Output[ENPY][ID2], DE_Status[ID2],
                         Gamma_m1, _Gamma_m1, CheckMinPres_No, NULL_REAL, DualEnergySwitch );
#     endif // #ifdef DUAL_ENERGY


//    check the negative density and energy
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CPU_CheckNegative(Output[DENS][ID2]) )
         Aux_Message( stderr, "WARNING : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      Output[DENS][ID2], __FILE__, __LINE__, __FUNCTION__ );

      if ( CPU_CheckNegative(Output[ENGY][ID2]) )
         Aux_Message( stderr, "WARNING : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      Output[ENGY][ID2], __FILE__, __LINE__, __FUNCTION__ );
#     endif

   } // i,j,k

} // FUNCTION : CPU_FullStepUpdate



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_StoreFlux
// Description :  Store the inter-patch fluxes for the AMR fix-up operation
//
// Parameter   :  Output   : Array to store the inter-patch fluxes
//                FC_Flux  : Array storing the face-centered fluxes
//                           --> Size is assumed to be N_FL_FLUX^3
//-------------------------------------------------------------------------------------------------------
void CPU_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL]  )
{

   int Face, ID1, ID2[9];

   for (int m=0; m<PS2; m++)
   for (int n=0; n<PS2; n++)
   {
      ID1    = m*PS2 + n;

      ID2[0] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + 0;
      ID2[1] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS1;
      ID2[2] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS2;
      ID2[3] = (  m*N_FL_FLUX +   0)*N_FL_FLUX + n;
      ID2[4] = (  m*N_FL_FLUX + PS1)*N_FL_FLUX + n;
      ID2[5] = (  m*N_FL_FLUX + PS2)*N_FL_FLUX + n;
      ID2[6] = (  0*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[7] = (PS1*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[8] = (PS2*N_FL_FLUX +   m)*N_FL_FLUX + n;

      for (int t=0; t<9; t++)
      {
         Face = t/3;

         for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Array[t][v][ID1] = FC_Flux[ ID2[t] ][Face][v];
      }
   }

} // FUNCTION : CPU_StoreFlux


#if ( COORDINATE == CYLINDRICAL )
// //-------------------------------------------------------------------------------------------------------
// Function    :  CurviFluxGrad
// Description :  get transverse flux gradient in curvilinear coordinate, used in FullStepUpdate 
//
// Parameter   :  Flux         : cell centered position at the corner cell of that patch, expect Corner_Array[3]
//                dF           : 
//                x_pos        : xyz position of the cell
//                face_pos     : xyz position of cell faces
//                v            : from 0 - (NCOMP_TOTAL-1)
//
// NOTE        :  could extend to all coordinate, but be careful of face_pos[][]
//-------------------------------------------------------------------------------------------------------
void CurviFluxGrad( real dF[][NCOMP_TOTAL], const real x_pos[] ) {

   // pre-calculate geo-factor
   const real _x1    = (real)1.0/x_pos[0];
   
   for (int d=0; d<3; d++)
   for (int v=0; v<NCOMP_TOTAL; v++) {
      if (v==MOMY) dF[d][v] *= SQR(_x1) ;
      else         dF[d][v] *= _x1      ;
   }
   
#  ifdef GAMER_DEBUG
   for (int d=0; d<3; d++)
   for (int v=0; v<NCOMP_TOTAL; v++) {
      if ( ! isfinite( dF[d][v] ) ) {
         Aux_Message( stderr, "WARNING : dF NaN'ed at (d, v) = (%d, %d) at file <%s>, line <%d>, function <%s>\n",
                      d, v, __FILE__, __LINE__, __FUNCTION__ );
      } // if ( !isfinite() ) 
   }
#  endif

}


// //-------------------------------------------------------------------------------------------------------
// Function    :  GetFullStepGeoSource
// Description :  get geometry source term that is used in full step update
//                use the predicted flux to calculate the corrected primitive var at half step
//                use this corrected half step pri var to find GeoSource for full step update
//
// Parameter   :  ConInput     : 
//                GeoSource    : 
//                dF           : 
//                x_pos        : 
//                dt_dh2       : 
//                dt_2         : dt/2
//
//-------------------------------------------------------------------------------------------------------
void GetFullStepGeoSource( const real ConInput[][ FLU_NXT*FLU_NXT*FLU_NXT ], real* GeoSource, 
                           const real dF[][NCOMP_TOTAL], const real* x_pos, const real* dt_dh2, const real dt_2, 
                           const real Gamma_m1, const real MinPres, const int ID3 ) {
                             
   real ConVar_Buffer[NCOMP_TOTAL], PriVar_Buffer[NCOMP_TOTAL];
   const bool NormPassive_No  = false; 
   const bool JeansMinPres_No = false;
   
   for (int v=0; v<NCOMP_TOTAL; v++) ConVar_Buffer[v] = ConInput[v][ID3] ;
   
   CPU_Con2Pri( ConVar_Buffer, PriVar_Buffer, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, 
                JeansMinPres_No, NULL_REAL );
   
   GeometrySourceTerm( PriVar_Buffer, x_pos, GeoSource );
   
   for (int v=0; v<NCOMP_TOTAL; v++) {
      ConVar_Buffer[v] -= dF[0][v]*dt_dh2[0] + dF[1][v]*dt_dh2[1] + dF[2][v]*dt_dh2[2] ;
      ConVar_Buffer[v] += GeoSource[v] * dt_2 ;
   }
   
   CPU_Con2Pri( ConVar_Buffer, PriVar_Buffer, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, 
                JeansMinPres_No, NULL_REAL );
   
   GeometrySourceTerm( PriVar_Buffer, x_pos, GeoSource );

#ifdef COOLING   
   
   real cool_rate;
   CoolingFunc(&cool_rate, PriVar_Buffer, x_pos);
   
   // ### note that GeoSource now includes both GeoSource and Cooling
   GeoSource[ENGY] -= cool_rate ; 
   
#endif
}

#endif // #if (COORDINATE == CYLINDRICAL)

#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

