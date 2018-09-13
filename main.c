//
//  main.c
//  HW3_Zhi_Su
//
//  Created by 苏治 on 2018/9/7.
//  Copyright © 2018年 治 苏. All rights reserved.
//

#include <petscvec.h>
#include <petscis.h>
#include <petscsys.h>
#include <petscviewer.h>

int main(int argc, const char * argv[]) {
    Vec            x,b,u;               /* vectors */
    PetscReal      norm;
    PetscErrorCode ierr;
    PetscScalar    one = 1.0,two = 2.0,three = 3.0,dots[3],dot;
    Ksp            ksp;
    Mat            A;
    PetscInt       n=100,its;
    PetscErrorCode ierr;
    PetscMPIInt    rank,size;
    
    
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscInitialize(&argc,&argv,(char*)0,help);
    PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);
    CHKERRQ(ierr);

    
    
    MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,&A);
    MatSetFromOptions(A);
    MatSetUp(A);
    /* (code to assemble matrix not shown) */
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE, n);
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
    VecDuplicate(x,&u);
    /*
    value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    for (i=1; i<n-1; i++) {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);
    }
    i= n - 1; col[0] = n - 2; col[1] = n - 1;
    MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
    MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    */
    /* (code to assemble RHS vector not shown)*/
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b, x);
    KSPDestroy(ksp);
    PetscFinalize();
    return 0;
}
