#include "GAMER.h"

// for POPIII
#include "TestProb.h"
#include <math.h>

#ifdef COOLING
extern double  Time[NLEVEL];
extern void    CoolingFunc(real cool_rate, const real PriVar[], const real x_pos[]);


//-------------------------------------------------------------------------------------------------------
// Function    :  CoolingFunc
// Description :  get cooling rate; 
//
// Parameter   :  cool_rate    : 
//                PriVar       : 
//
// NOTE        :  
//-------------------------------------------------------------------------------------------------------
void CoolingFunc(real cool_rate, const real PriVar[], const real x_pos[]) {
   
   // Boltzmann R in the unit of popIII setting
   const double t_orbit     = 0.79 ;                  // outer orbital time
   const double R           = 4.64952804093003e+0 ; 
   
   const double t_curr      = Time[0];
   const double t_relax     = 5*t_orbit ; 
   
   
   const double GM          = ExtAcc_AuxArray[3] ;
   const double rho         = PriVar[DENS]; 
   const double pres        = PriVar[ENGY];
   const double ie          = PriVar[4]/(GAMMA-1.0) ;
   const double T           = pres / ( rho * R ) ;
   const double v_abs       = SQRT( SQR(PriVar[MOMX]) + SQR(PriVar[MOMY]) + SQR(PriVar[MOMZ]) ) ;
   const double sph_rad     = SQRT( SQR(x_pos[0]) + SQR(x_pos[2]) );
   
   //const double tau_dyn     = SQRT( CUBE(sph_rad) / GM );
   const double tau_dyn     = sph_rad / v_abs ;
   
   real cool_time = 1.0 * tau_dyn;
   
   // check if temperature is smaller than 100K
   if (T>= 100)   cool_rate = ie / cool_time ;
   else           cool_rate = TINY_NUMBER ;
   
   if (t_curr < t_relax) cool_rate *= FABS( 1.0 - (t_relax-t_curr)/t_relax ) ; 
   

   
} // FUNCTION: CoolingFunc


#endif // COOLING
