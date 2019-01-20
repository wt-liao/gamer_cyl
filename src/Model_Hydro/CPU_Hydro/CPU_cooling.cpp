#include "GAMER.h"

// for POPIII
#include "TestProb.h"
#include <math.h>

#ifdef COOLING
extern void CoolingFunc(real cool_rate, const real PriVar[], const real x_pos[]);


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
   
   // This script assume GM=1; const_R=1 
   const double GM          = 1.0;
   const double const_R     = 1.0; 
   const double T_0         = 1.0;
   const double R_0         = 1.0;
   
   const double tau_dyn     = 1.0; // <= SQRT(GM/R_0^3)
   
   const double rho         = PriVar[DENS];
   const double T           = PriVar[ENGY]/rho/const_R ;
   
   const double beta        = 0.5;
   const double cool_time   = beta * tau_dyn;
   
   cool_rate = rho/(GAMMA-1)*(T-T_0)/cool_time;
   cool_rate = FMAX(cool_rate, 0.0);
   
} // FUNCTION: CoolingFunc




#endif // COOLING