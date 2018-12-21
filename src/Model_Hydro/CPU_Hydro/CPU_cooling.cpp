#include "GAMER.h"

// for POPIII
#include "TestProb.h"
#include <math.h>

#ifdef COOLING
extern void CoolingFunc(real cool_rate, const real PriVar[], const real x_pos[]);

// for POPIII
static double rate_k4_func(const real T);
static double rate_k5_func(const real T);
static double f_HII_func(const real n_cgs, const real T);
static double f_H2_func(const real n_cgs, const real T, const real tau_dyn_cgs);
static double H2_cool_func(const real rho_cgs, const real T, const real f_H2, const real m_H_cgs, const double X);
static double CIE_cool_func(const real rho_cgs, const real T, const real f_H2, const double X);
static double Ly_cool_func(const real T, const real n_e, const double n_HI);
static double Brem_cool_func(const real T, const real n_e, const double n_HII);

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
   
   // this is an example for const cooling time
   /*
   real cool_time = 0.002;
   
   real ie;
   ie = PriVar[4]/(GAMMA-1.0) ; 
   
   cool_rate = ie / cool_time ;
   */
   
   // POPIII
   const double m_H         = 8.41743046550226e-52; // dimensionless
   const double R           = 4.64952804093003e+0 ; 
   const double GAMMA_1     = GAMMA - (real)1.0 ;
   
   const double rho_unit    = 5.93924146044389e-13;
   const double mass_unit   = 1.98841586000000e+27;
   const double L_unit      = 1.49597870750767e+13;
   const double time_unit   = 5.02280842159966e+09;
   const double rate_unit   = mass_unit/L_unit/CUBE(time_unit); //erg cm^{-3} s^(-1) = g*cm2*s-2*cm-3*s-1 = g*cm-1*s-3
      
   const double GM          = ExtAcc_AuxArray[3] ;
   const double sph_rad     = SQRT( SQR(x_pos[0]) + SQR(x_pos[2]) );
   const double tau_dyn     = SQRT( CUBE(sph_rad) / GM );
   
   const double rho         = PriVar[DENS]; 
   const double pres        = PriVar[ENGY];
   const double T           = pres / ( rho * R ) ;
   const double n           = rho/(2.0*m_H) ;
   
   // convert to cgs
   const double rho_cgs     = rho*rho_unit ;
   const double m_H_cgs     = m_H * mass_unit;
   const double n_cgs       = n / CUBE(L_unit);
   const double tau_dyn_cgs = tau_dyn * time_unit; 
   
   //###
   const double X           = 1.0;
   
   //double f_HI, f_HII, f_H2, n_HII, n_e, n_HI;
   double H2_cool, CIE_cool, Ly_cool, brem_cool; // cooling rate
   
   // composition
   const double f_HII  = f_HII_func(n_cgs, T);
   const double f_H2   = f_H2_func(n_cgs, T, tau_dyn_cgs) ;
   const double f_HI   = 1.0 - f_HII - f_H2;
   
   const double n_HII  = n_cgs * f_HII ;
   const double n_e    = n_HII;
   const double n_HI   = n_cgs * f_HI;
   
   // cooling rate
   H2_cool   = H2_cool_func(rho_cgs, T, f_H2, m_H_cgs, X);
   CIE_cool  = CIE_cool_func(rho_cgs, T, f_H2, X);
   Ly_cool   = Ly_cool_func(T, n_e, n_HI);
   brem_cool = Brem_cool_func(T, n_e, n_HII) ;
   
   //### reduce the cool_rate by 0.2 in first outer orbit
   //### make sure this commented out after 3 orbit
   H2_cool  *= 0.2;
   CIE_cool *= 0.2;

   
   // dimensionless cool rate
   cool_rate = (H2_cool + CIE_cool + Ly_cool + brem_cool) / rate_unit ;
   
   // for debug
   /*
   double ie = PriVar[4]/(GAMMA-1.0) ;
   double cool_time = ie / cool_rate ;
   
   if ( n_cgs>1e13 && T<4000 )
      Aux_Message(stdout, "At (rho, n, T) = (%5.2e, %5.2e, %6.0f). Cooling time = %5.2e. (H2, CIE, Ly, Brem) = (%4.1e, %4.1e, %4.1e, %4.1e). \n", 
                  rho_cgs, n_cgs, T, cool_time/tau_dyn, H2_cool, CIE_cool, Ly_cool, brem_cool);
   */

   
} // FUNCTION: CoolingFunc


// POPIII
//-------------------------------------------------------------------------------------------------------
// Function    :  rate_k4 in cgs
//-------------------------------------------------------------------------------------------------------
double rate_k4_func(const real T) {
   double rate_k4 = 5.5e-29 * POW(T, -1) ;
   
   return rate_k4 ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  rate_k5 in cgs
//-------------------------------------------------------------------------------------------------------
double rate_k5_func(const real T) {
   double rate_k5 = 6.5e-7 * POW(T, -0.5) * EXP(-52000/T) * (1-EXP(-6000/T))  ;   
   
   return rate_k5 ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  f_HII 
//-------------------------------------------------------------------------------------------------------
double f_HII_func(const real n_cgs, const real T) {
   const double consts  = 2.86457498307206e+09 ;     // m_e*k_B/h^2 in units of K^(-1) cm^(-2)
   const double T_ratio = 1.57821462222688e+05 / T ;
   
   double c1 = POW(2*M_PI*consts*T, 1.5) / n_cgs * EXP(- T_ratio); 
   double x  = 0.5*( -c1 + SQRT(c1*c1 + 4.0*c1) );
      
   return x ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  f_H2
//-------------------------------------------------------------------------------------------------------
double f_H2_func(const real n_cgs, const real T, const real tau_dyn_cgs) {
   const double k4    = rate_k4_func(T);
   const double k5    = rate_k5_func(T); 
   const double coeff = k5/(2*k4);
   
   double n_HI, n_H2_estimate, n_H2, tau_chem;
   
   n_HI          = 0.5*( -coeff + SQRT(coeff*coeff + 4*coeff*n_cgs) ) ;
   n_H2_estimate = (k4/k5)*SQR(n_HI) ;
   
   tau_chem = 1 / ( k4*SQR(n_cgs) ) ;
   n_H2     = n_H2_estimate * EXP(-tau_chem / tau_dyn_cgs);
   
   return FMAX(2*n_H2 / n_cgs, 1.0e-8) ;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  H2_cool rate in erg cm^(-3) s^(-1)
//-------------------------------------------------------------------------------------------------------
double H2_cool_func(const real rho_cgs, const real T, const real f_H2, const real m_H_cgs, const double X) {
   const double T3 = T*1e-3;
   double rate;
   
   // rate in erg g-1 s-1
   rate = X*f_H2/m_H_cgs*( (9.5e-22*POW(T3,3.76))/(1+0.12*POW(T3,2.1)) * EXP(-POW(0.13/T3,3)) + 3e-24*EXP(-0.51/T3)
                          + 6.7e-19*EXP(-5.86/T3) + 1.6e-18*EXP(-11.7/T3) ) ;    
   // rate = rate * rho
   rate *= rho_cgs ;
   
   return rate;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  CIE_cool rate in erg cm^(-3) s^(-1)
//-------------------------------------------------------------------------------------------------------
double CIE_cool_func(const real rho_cgs, const real T, const real f_H2, const double X) {
   
   // rate in erg g-1 s-1
   double rate = 7.2e-2*rho_cgs*POW(T, 4)*X*f_H2 ;
   // rate = rate * rho
   rate *= rho_cgs;
   
   return rate;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  Ly_cool rate in erg cm^(-3) s^(-1)
//-------------------------------------------------------------------------------------------------------
double Ly_cool_func(const real T, const real n_e, const double n_HI) {
   const double T5 = T*1e-5;
   
   // rate in erg cm-3 s-1
   double rate = 7.5e-19/(1+SQRT(T5)) * EXP(-118348/T) *n_e * n_HI ;
   
   return rate;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  Brem_cool rate in erg cm^(-3) s^(-1)
//-------------------------------------------------------------------------------------------------------
double Brem_cool_func(const real T, const real n_e, const double n_HII) {
   
   // rate in erg cm-3 s-1
   double rate = 1.43e-27*SQRT(T) * (1.1 + 0.34*EXP(- SQR(5.5-log10(T))/3.0 ) ) * n_e * n_HII ;
   
   return rate;
}


#endif // COOLING
