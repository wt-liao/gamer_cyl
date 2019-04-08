#include "GAMER.h"


//### decalre Grackle_Load_Alpha_Table

#if (defined SUPPORT_GRACKLE) && (defined GRACKLE_H2_SOBOLEV) 
//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Load_Alpha_Table
// Description :  Load table of H2 abosorption coefficient
//-------------------------------------------------------------------------------------------------------
void Grackle_Load_Alpha_Table(){
   
   // read in 
   FILE * T_Table;
   FILE * Alpha_Table; 
   
   //### T_Table should be in Ln(T) space
   T_Table     = fopen( "H2_Op_T" ,    "rb" );
   Alpha_Table = fopen( "H2_Op_Alpha", "rb" );
   
   if (T_Table == NULL && MPI_Rank == 0)
      Aux_Message(stderr, "ERROR: cannot find <H2_Op_T>: H2 Opacity T-table . \n");
   if (Alpha_Table == NULL && MPI_Rank == 0)
      Aux_Message(stderr, "ERROR: cannot find <H2_Op_Alpha>: H2 Opacity Alpha-table. \n");
   
   fseek (T_Table ,     0 , SEEK_END);
   fseek (Alpha_Table , 0 , SEEK_END);
   
   long Table_T_Size     = ftell(T_Table)    ;
   long Table_Alpha_Size = ftell(Alpha_Table);
   
   if (Table_T_Size != Table_Alpha_Size)
      Aux_Message(stderr, "Error in H2 Opacity: size of T-table=%d and Alpha-table=%d. \n ", 
                  Table_T_Size, Table_Alpha_Size) ;
   
   const int N_elem = int( Table_T_Size/sizeof(double) ); 
   
   if (MPI_Rank == 0)
      Aux_Message(stdout, "H2 Opacity NOTE: %d elements in H2 Opacity Table. \n", N_elem); 
   
   // allocate memory for T-table and Alpha-table
   //### declare H2_Op_T_Tabl, H2_Op_Alpha_Table in global.h -> Done!
   H2_Op_T_Table     = new double[N_elem]; 
   H2_Op_Alpha_Table = new double[N_elem];
   
   // save Table to array
   rewind(T_Table);
   rewind(Alpha_Table);
   
   fread(H2_Op_T_Table,     sizeof(double), N_elem, T_Table    ); 
   fread(H2_Op_Alpha_Table, sizeof(double), N_elem, Alpha_Table); 
   
   //### rescale H2_Op_T_Table /= T_unit 
   //### -> currently make no difference since T_unit = 1
   
   // close file
   fclose(T_Table); 
   fclose(Alpha_Table);
   
   //### allocate Grackle_T_Start, Grackle_dT in global.h. Done!
   Grackle_T_Start = H2_Op_T_Table[0]; 
   Grackle_T_End   = H2_Op_T_Table[N_elem -1] ;
   Grackle_dT      = H2_Op_T_Table[1] - H2_Op_T_Table[0] ;
   
   //### free H2_Op_T_Table, H2_Op_Alpha_Table in Grackle_End(). Done!
   
}

#endif