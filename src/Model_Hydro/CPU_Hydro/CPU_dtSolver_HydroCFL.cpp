#include "GAMER.h"

#if ( !defined GPU  &&  MODEL == HYDRO )

#ifdef COOLING
extern void CoolingFunc(real* cool_rate, const real PriVar[], const real x_pos[]);
#endif


//-----------------------------------------------------------------------------------------
// Function    :  CPU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) required for the hydro solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//
// Parameter   :  dt_Array     : Array to store the minimum dt in each target patch
//                Flu_Array    : Array storing the prepared fluid data of each target patch
//                Corner_Array : Array storing the physical corner coordinates of each patch
//                NPG          : Number of target patch groups
//                dh           : Grid size
//                Safety       : dt safety factor
//                Gamma        : Ratio of specific heats
//                MinPres      : Minimum allowed pressure
//
// Return      :  dt_Array
//-----------------------------------------------------------------------------------------
void CPU_dtSolver_HydroCFL( real dt_Array[], const real Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const double Corner_Array[][3],
                            const int NPG, const real dh[], const real Safety, const real Gamma, const real MinPres )
{

   const bool CheckMinPres_Yes = true;
   const int  NPatch           = 8*NPG;
   const real Gamma_m1         = Gamma - (real)1.0;
#  ifdef COOLING
   const real _safety_cool = 2.0;
#  endif

   real fluid[NCOMP_FLUID], _Rho, Vx, Vy, Vz, Pres, Cs, CurrCFL, MaxCFL;

// loop over all patches
#  pragma omp parallel for private( fluid, _Rho, Vx, Vy, Vz, Pres, Cs, CurrCFL, MaxCFL ) schedule( runtime )
   for (int p=0; p<NPatch; p++)
   {
      MaxCFL = (real)0.0;
      real _dh[3] = { (real)1.0/dh[0], (real)1.0/dh[1], (real)1.0/dh[2] };
      real x_pos[3] ;
      int ID; 
#     ifdef COOLING
      real _dt_cool, cool_rate, PriVar[NCOMP_TOTAL];
#     endif

      for (int k=0; k<PS1; k++)
      for (int j=0; j<PS1; j++)
      for (int i=0; i<PS1; i++) 
      {
         ID = (k*PS1 +j)*PS1 +i;
         for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = Flu_Array[p][v][ID];
         
#        if ( COORDINATE == CYLINDRICAL )
         x_pos[0] = Corner_Array[p][0] + i*dh[0] ;
         x_pos[1] = Corner_Array[p][1] + j*dh[1] ;
         x_pos[2] = Corner_Array[p][2] + k*dh[2] ;
         
         const real radius = Corner_Array[p][0] + i*dh[0] ;
         _dh[1] = (real)1.0/ ( dh[1]*radius ) ;
#        endif

        _Rho  = (real)1.0 / fluid[DENS];
         Vx   = FABS( fluid[MOMX] )*_Rho;
         Vy   = FABS( fluid[MOMY] )*_Rho;
         Vz   = FABS( fluid[MOMZ] )*_Rho;
         Pres = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                 Gamma_m1, CheckMinPres_Yes, MinPres );
         Cs   = SQRT( Gamma*Pres*_Rho );

#        if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU  ||  FLU_SCHEME == WAF )
         CurrCFL = FMAX( (Vx+Cs)*_dh[0], (Vy+Cs)*_dh[1] );
         CurrCFL = FMAX( (Vz+Cs)*_dh[2], CurrCFL );
         MaxCFL  = FMAX( CurrCFL, MaxCFL );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
         CurrCFL = (Vx+Cs)*_dh[0] + (Vy+Cs)*_dh[1] + (Vz+Cs)*_dh[2];
         MaxCFL  = FMAX( CurrCFL, MaxCFL );
         
         //### may also need cool_dt_safty: _dt_cool * safty (safty > 1; probably 10)
#        ifdef COOLING
         PriVar[0] = fluid[DENS];
         PriVar[1] = Vx;
         PriVar[2] = Vy;
         PriVar[3] = Vz;
         PriVar[4] = Pres;
         
         CoolingFunc(&cool_rate, PriVar, x_pos);
         _dt_cool = cool_rate * (Gamma-1.0) / Pres ;
         MaxCFL   = FMAX(_safety_cool*_dt_cool, MaxCFL) ;
#        endif // #ifdef COOLING
         
#        endif
         
      } // for (int t=0; t<CUBE(PS1); t++)

      dt_Array[p] = Safety/MaxCFL;

   } // for (int p=0; p<NPatch; p++)

} // FUNCTION : CPU_dtSolver_HydroCFL



#endif // #if ( !defined GPU  &&  MODEL == HYDRO )
