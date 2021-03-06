#include <petscksp.h>
int main(int argc, char **args) {
    Vec            x,b,u;
    PetscReal      norm,tol=1000.*PETSC_MACHINE_EPSILON;
    PetscScalar    one=1.0,negone=-1.0,matvalue[3];
    KSP            ksp;
    PC             pc;
    Mat            A;
    PetscInt       i,n=100,colum[3],its;
    PetscMPIInt    rank,size;
    
    PetscInitialize(&argc,&args,(char *)0,NULL);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
    MatSetFromOptions(A);
    MatSetUp(A);
    //built matrix A for sparse form//
        for (i=0; i<n; i++) {
        if (i==0) {
            colum[0]=0;
            colum[1]=1;
            matvalue[0]=2.0;
            matvalue[1]=-1.0;
            MatSetValues(A,1,&i,2,colum,matvalue,INSERT_VALUES);
        }
        else if (i>0&&i<n-1){
            matvalue[0]=-1.0;
            matvalue[1]=2.0;
            matvalue[2]=-1.0;
            colum[0]=i-1;
            colum[1]=i;
            colum[2]=i+1;
            MatSetValues(A,1,&i,3,colum,matvalue,INSERT_VALUES);
        }
        else if (i==n-1){
            colum[0]=i-1;
            colum[1]=i;
            matvalue[0]=-1.0;
            matvalue[1]=2.0;
            matvalue[2]=-1.0;
            MatSetValues(A,1,&i,2,colum,matvalue,INSERT_VALUES);
        }
    }
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    //build RHS B//
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,n);
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
    VecDuplicate(x,&u);
    VecSet(b,one);
    //KSP
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCJACOBI);
    KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp,b,x);
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    //check ||AX-b||
    MatMult(A,x,u);
    VecAXPY(x,negone,u);
    VecNorm(x,NORM_2,&norm);
    KSPGetIterationNumber(ksp,&its);
    if(norm>tol){
    PetscPrintf(PETSC_COMM_WORLD,"residual is %g, Iterations %D\n",norm,its);
   }

    VecDestroy(&x);
    VecDestroy(&b);
    VecDestroy(&u);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    PetscFinalize();
    return 0;
}

