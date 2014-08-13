/*! @file VecMarket_to_PETScBin.c
 *        @brief : This routine converts Vector in matrix market file to PETSc Binary format
*/

/*
* For this implementation, we have used variours PETSc routines *
* For more information, check http://www.mcs.anl.gov/petsc/petsc-as/ and specifically  ksp/ex72.c & mat/ex72.c 
*/

/* PETSc Headers */
#include "petsc.h"
#include "petscmat.h"

int main(int argc,char **args)
{

  /*Petsc Vector Object*/
  Vec 	      b;
  /*Input and Output file names*/
  char        filein[128],fileout[128],buf[128];

  int         i,m,n,ierr;

  PetscScalar val;
  /*File Handle for Reading matrix market file */
  FILE*       file;

  /*Petsc viewer for reading and writing binary vector */
  PetscViewer view;
  PetscLogDouble t1,t2,elapsed_time;

  /*Initialize PETSc library */
  PetscInitialize(&argc,&args,(char *)0,PETSC_NULL);
  /*number of comments line in matrix market file*/

  int numComments=0;

  ierr = PetscGetTime(&t1); CHKERRQ(ierr);

  /* Open Matrix Market File for Vector */
  ierr = PetscOptionsGetString(PETSC_NULL,"-fin",filein,127,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscFOpen(PETSC_COMM_SELF,filein,"r",&file);CHKERRQ(ierr);

/* Just count the comment lines in the file */
  while(1)
  {
  	fgets(buf,128,file);
        /*If line starts with %, its a comment */
        if(buf[0] == '%')
	{
	   printf("\n IGNORING COMMENT LINE : IGNORING....");
	   numComments++; 
	}
	else
	{
	   /*Set Pointer to Start of File */
	   fseek(file, 0, SEEK_SET );
           int num = numComments;

	   /* and just move pointer to the entry in the file which indicates row nums, col nums and non zero elements */
	   while(num--)
	   	fgets(buf,128,file);
	   break;
	}
  }

  /* Reads Size of Vector form Matrix Market file */
  fscanf(file,"%d %d\n",&m,&n);
  printf ("M = %d, N = %d\n",m,n);

  /*Create PETSc Vector Objects */
  ierr = VecCreate(PETSC_COMM_WORLD,&b);
  ierr = VecSetSizes(b,m,PETSC_DECIDE);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
  
  /*Now Reads Vector From File and Set values to PETSc Object */
  for (i=0; i<m; i++) {
	    fscanf(file,"%le \n",&val);
    	ierr = VecSetValues(b,1, &i, &val, INSERT_VALUES);
  }

  /*Assmeble the Vector */
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  fclose(file);

  /* Create PETScViewer for Writing Binary Vector */
  ierr = PetscOptionsGetString(PETSC_NULL,"-fout",fileout,127,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileout,FILE_MODE_WRITE,&view);CHKERRQ(ierr);
  /*Dump the Vector to binary format*/
  ierr = VecView(b,view);CHKERRQ(ierr);
  
  /*Destroy the Data structure */
  ierr = PetscViewerDestroy(&view);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);

  ierr = PetscGetTime(&t2);CHKERRQ(ierr);
  elapsed_time = t2 - t1;     
  ierr = PetscPrintf(PETSC_COMM_SELF,"ELAPSED TIME %g\n",elapsed_time);CHKERRQ(ierr);

  /*End Program*/
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
