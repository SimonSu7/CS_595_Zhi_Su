static char help[] = "Builds a parallel vector with 1 component on the first processor, 2 on the second, etc.\n\ Then each processor adds one to all elements except the last rank.\n\n";

#include <petscsys.h>
#include <petscksp.h>

int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       i,j,N;
  PetscScalar    one = 1.0;
  Vec            x;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);
  if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,rank+1,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetSize(x,&N);CHKERRQ(ierr);
  ierr = VecSet(x,one);CHKERRQ(ierr);

  for (i=0,j=0; i<N-rank; i++,j++) {
    ierr = VecSetValues(x,1,&i,&one,INSERT_VALUES);CHKERRQ(ierr);
    if(j==one-1){
       j=-1;
       one++;
    }
  }


  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
        

  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
    
  ierr = PetscFinalize();
  return ierr;
} 
