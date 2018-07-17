#include "GAMER.h"

#if (COORDINATE == CYLINDRICAL)
#ifdef GRAVITY

extern rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;    


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_CylKernel
// Description :  Calculate Kernel function in k-space
//-------------------------------------------------------------------------------------------------------
void Init_CylKernel(){
   
   const double *dh      = amr->dh[0]; 
   const double dh_cube  = dh[0]*dh[1]*dh[2] ;
   const int FFT_Size[3] = { NX0_TOT[0], NX0_TOT[1], NX0_TOT[2]*2 }; // FFT_Size[0] is redundunt

   double x, xp, y, z ; // (x, y, z) <-> (r, phi, z)
   int ID_planX, ID_planYZ, kk, ii, iip ;
      
# ifdef SERIAL
   const int local_nx        = NX0_TOT[0];
   const int local_nxp       = NX0_TOT[0]; 
   const int local_ny        = 2*(FFT_Size[1]/2+1);
   const int local_nz        = FFT_Size[2]; 
   const int local_nx_start  = 0;    
   const int local_nxp_start = 0 ;
   //const int RANK_I_TOT      = 1;
   //const int RANK_IP_TOT     = 1;
# else
   //### arrange the cores to be MPI_NRANK * 1
   //### this arrangement will have memory cap by the size of 3D array
   //### NEED a scheme to determine RANK_I_TOT, RANK_IP_TOT for general case 
   
   const int RANK_I_TOT  = MPI_NRank;
   const int RANK_IP_TOT = 1;       // CEIL ( (2*NX0_TOT[0]*NX0_TOT[1]*NX0_TOT[2])/memory_cap );
   
   // MPI_RANK = RANK_IP*(RANK_I_TOT) + RANK_I 
   const int RANK_IP         = int(MPI_Rank/RANK_I_TOT);          // const int RANK_IP = 0;
   const int RANK_I          = MPI_Rank % RANK_I_TOT ;            // const int RANK_I  = MPI_Rank;    
   const int global_nx_unit  = ceil(NX0_TOT[0]/RANK_I_TOT );
   const int global_nxp_unit = ceil(NX0_TOT[0]/RANK_IP_TOT);
   const int local_nx_unit   = ceil(global_nx_unit / RANK_IP_TOT);
   const int local_nxp_unit  = ceil(global_nxp_unit/ RANK_I_TOT );
   
   int global_nx, global_nxp, local_nx, local_nxp ;
   if (RANK_I  != RANK_I_TOT-1)  global_nx  = global_nx_unit;   
   else                          global_nx  = NX0_TOT[0] - global_nx_unit *(RANK_I_TOT -1);
   if (RANK_IP != RANK_IP_TOT-1) global_nxp = global_nxp_unit;
   else                          global_nxp = NX0_TOT[0] - global_nxp_unit*(RANK_IP_TOT-1);
   if (RANK_I  != RANK_I_TOT-1)  local_nxp  = local_nxp_unit;   
   else                          local_nxp  = global_nxp_unit - local_nxp_unit*(RANK_I_TOT -1);
   if (RANK_IP != RANK_IP_TOT-1) local_nx   = local_nx_unit;
   else                          local_nx   = global_nx_unit  - local_nx_unit *(RANK_IP_TOT-1);

   const int  global_nx_start  = RANK_I *global_nx_unit ;
   const int  global_nxp_start = RANK_IP*global_nxp_unit;
   const int  local_ny         = 2*(FFT_Size[1]/2+1);
   const int  local_nz         = FFT_Size[2]; 
   const long slab_size        = local_ny * local_nz ;
   
   // init memory for MPI
   Init_MemAllocate_CylPoisson(local_nx, local_nxp, global_nx, global_nxp, slab_size);
   
# endif // ifdef SERIAL, else...
   
   // new MPI Comm
   MPI_Comm_split(MPI_COMM_WORLD, RANK_I,  MPI_Rank, &rank_i_comm );
   MPI_Comm_split(MPI_COMM_WORLD, RANK_IP, MPI_Rank, &rank_ip_comm);
   
   if (rank_i_comm == MPI_COMM_NULL) Aux_Error(ERROR_INFO, "new MPI_Comm initiation failed... \n") ;

   
   // 1.0 build up kernel for FFT
   // 1.1 allocate Kernel K - 
   //   a 1D array (a flattened 2D array), each component is a pointer to another 1D array (a flattened 2D array)
   Aux_AllocateArray2D(KernelFuncK, global_nx*global_nxp, slab_size) ;
   
   // 1.2 build up Kernel in real space
   for (int i=0;  i <global_nx;  i++)  { ii  = i+ global_nx_start;  x  = amr->BoxEdgeL[0] + (ii +0.5)*dh[0];
   for (int ip=0; ip<global_nxp; ip++) { iip = ip+global_nxp_start; xp = amr->BoxEdgeL[0] + (iip+0.5)*dh[0];
                                         ID_planX = i*global_nxp + ip;
   for (int k=0;  k <local_nz;  k++)   { z = ( k <= NX0_TOT[2] ) ? k*dh[2] : (FFT_Size[2]-k)*dh[2] ;
   for (int j=0;  j <local_ny;  j++)   { y = j*dh[1];
      
      ID_planYZ = k*local_ny + j;
      real denominator = SQRT( SQR(x-xp) + (real)2.0*x*xp*((real)1.0 - COS(y)) + SQR(z) ) ;
      
      if (denominator != 0.0) 
         KernelFuncK[ID_planX][ID_planYZ] = (real) -1.0*dh_cube / denominator; 
      else  
         KernelFuncK[ID_planX][ID_planYZ] = (real) 0.0 ;    // mesh does not see itself
         
   }}}}
   
   
   // 2.0 FFT Kernel for each (i, ip)
   for (int i=0;  i<global_nx;   i++){
   for (int ip=0; ip<global_nxp; ip++){
      
      ID_planX = i*global_nxp + ip;
      rfftwnd_one_real_to_complex( FFTW_Plan, KernelFuncK[ID_planX], NULL );      
   }}
      
} // Init_CylKernel


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_CylPoisson
// Description :  allocate memory needed for CylPoisson (promary design for MPI task)
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_CylPoisson(const int local_nx, const int local_nxp, const int global_nx, 
                                 const int global_nxp, const long slab_size){
   // need parameters: local_nx, local_nxp, slab_size
   
   const int  NSlab           = amr->NPatchComma[0][1]*PS1;   // number of slabs in each node 
   const int  PSSize          = PS1 * PS1;
   const long slab_size_hf    = slab_size / 2 ;
   const long local_nxp_total = local_nxp * NX0_TOT[1] * NX0_TOT[2];
   const long local_nx_total  = local_nx  * NX0_TOT[1] * NX0_TOT[2];
   const long local_nxp_slab  = local_nxp_total / PSSize;
   const long local_nx_slab   = local_nx_total  / PSSize;
   
   
   // 1.0 memory in CPU_CylPoissonSolver
   Aux_AllocateArray2D(RhoK, local_nxp, slab_size) ;
   Aux_AllocateArray2D(PhiK, local_nx , slab_size) ;
   
   // 2.0 memory in Patch2Slab
   SendBuf_Rho      = new real [ NSlab*PSSize ]; 
   SendBuf_IDPlanXp = new int  [ NSlab  ];
   SendBuf_IDPlanYZ = new long [ NSlab  ];
   RecvBuf_Rho      = new real [ local_nxp_total ]; 
   RecvBuf_IDPlanXp = new int  [ local_nxp_slab  ];
   RecvBuf_IDPlanYZ = new long [ local_nxp_slab  ];
   
   // 3.0 memory in Pot_Isolated
   SendBuf_RhoK_re  = new real [ local_nxp  * slab_size ] ;
   SendBuf_RhoK_im  = new real [ local_nxp  * slab_size ] ;
   RhoK_All_re      = new real [ global_nxp * slab_size ];
   RhoK_All_im      = new real [ global_nxp * slab_size ];
   PhiK_All_re      = new real [ global_nx  * slab_size_hf ] ; 
   PhiK_All_im      = new real [ global_nx  * slab_size_hf ] ;
   PhiK_local_re    = new real [ local_nx   * slab_size_hf ] ;
   PhiK_local_im    = new real [ local_nx   * slab_size_hf ] ;
   
   // 4.0 in Slab2Patch   
   SendBuf_Phi      = new real [ local_nx_total ]; 
   SendBuf_PID      = new long [ local_nx_slab  ];
   SendBuf_I        = new int  [ local_nx_slab  ];
   // can reuse the one in step 1.0
   RecvBuf_Phi      = new real [ NSlab*PSSize ] ;
   RecvBuf_PID      = new long [ NSlab        ] ;
   RecvBuf_I        = new int  [ NSlab        ] ;
   
}  //FUNCTION: Init_MemAllocate_CylPoisson

#endif // GRAVITY
#endif // COORDINATE == CYLINDRICAL