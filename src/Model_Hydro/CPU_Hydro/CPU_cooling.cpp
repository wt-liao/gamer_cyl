#include "GAMER.h"

extern void CoolingFunc(real cool_rate, const real PriVar);


//-------------------------------------------------------------------------------------------------------
// Function    :  CoolingFunc
// Description :  get cooling rate; 
//
// Parameter   :  cool_rate    : 
//                PriVar       : 
//
// NOTE        :  
//-------------------------------------------------------------------------------------------------------
void CoolingFunc(real cool_rate, const real PriVar) {
   
   // this is an example for const cooling time
   real cool_time = 10.0;
   
   real ie;
   ie = PriVar[4]/(GAMMA-1.0) ; 
   
   cool_rate = ie / cool_time ;
   
} // FUNCTION: CoolingFunc