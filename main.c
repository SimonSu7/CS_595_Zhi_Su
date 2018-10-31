//
//  main.c
//  CS 595 Final Project
//
//  Created by 苏治 on 2018/10/11.
//  Copyright © 2018年 治 苏. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define F "/Users/zhisu/Desktop/CS 595/pro2/1.csv"
#define OUT "/Users/zhisu/Desktop/CS 595/pro2/out.csv"
#define N 570
// X 就是取outer boundary 的判断个数
#define X 5
// M 就是line 被多少pixels 分割
#define M 1
// B is how many pixel in a line
#define B 50
//L level
#define L 121
void read_data();
void smooth();
void binary();
void boundary();
void subsample();
void chaincode();
void writeExcel();
void outer_boundary();

struct startpoint{
    int firsti;
    int firstj;
    int selecti;
    int selectj;
    int chaincode[N];
    int numofsub;
}start;
struct image{
    int data[N];
    int smoothed_data[N];
    int binary[N];
    int boundary[N];
    int subsample[N];
}guy[N];

int main(int argc, const char * argv[]) {
    read_data();
    smooth();
    binary();
    boundary();
    outer_boundary();
    subsample();
    writeExcel();
    chaincode();
}

void read_data()
{
    FILE *fp;
    int i,j = 0;
    fp=fopen(F,"r");
    //printf("\t");
    fseek(fp, 3L, SEEK_SET);   // 从文件第二行开始读取
    for(i = 0 ;i <N ; i++){
        for(j = 0 ;j < N ; j++)
        {
            fscanf(fp,"%d",&guy[i].data[j]);
            fseek(fp, 1L, SEEK_CUR);   /*fp指针从当前位置向后移动*/
            //printf("%d\t",guy[i].data[j]);
        }
        //puts("");
    }
   
    //puts("");
    //printf("%d",guy[6].data[9]);
    //puts("");
    fclose(fp);
}
// converlutional computation
void smooth(){
    int i,j,p,q;
    // m is average filter, n is part of the image(9*9); x is extended data h is double inversed 9*9 image matrix
    float r=0;
    puts("");
    //build new image matrix with 0 edge
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            guy[i].smoothed_data[j]=0;
        }
    }
    int x[N+8][N+8]={0};
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            x[i+4][j+4]=guy[i].data[j];
        }
    }
    
    /*show x matrix
    for (i=0; i<N+8; i++) {
        for (j=0; j<N+8;j++) {
            printf("%d ",x[i][j]);
        }
        puts("");
    }*/
    
    //build matrix filter
    /*for (i=0; i<9; i++) {
        for (j=0; j<9; j++) {
            m[i][j]=1/81;
        }
    }*/
    
    //apply numerical computing
    //1 step seperate the 9*9 matrix from image matrix(no use)
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            for (p=0; p<9; p++) {
                for (q=0; q<9; q++) {
                    r=r+x[i+p][j+q];
                }
            }
            guy[i].smoothed_data[j]=r/81;
            r=0;
            //printf("%d\t",guy[i].smoothed_data[j]);
        }
        //puts("");
    }
}
void binary(){
    int i,j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (guy[i].smoothed_data[j]>L) {
                guy[i].binary[j]=1;
            }
            else guy[i].binary[j]=0;
            //printf("%d ",guy[i].binary[j]);
        }
        //puts("");
    }
}
// second way to realize boundary
/*void boundary(){
    int i,j,bull = 0,newi=0,newj=0,n,x=0,y=0,z=0,h=0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            //printf("%d ",guy[i].boundary[j]);
            if (guy[i].binary[j]==1) {
                bull=1;
                start.firsti=i;
                start.firstj=j;
                break;
            }
        }
        if (bull==1) {
            break;
        }
    }
    i=start.firsti;
    j=start.firstj;
    while (i!=newi&&j!=newj) {
        x=0;
        y=0;
        z=0;
        h=0;
        for (n=0; n<X; n++) {
            x=x+guy[i].binary[j-n-1];
            y=y+guy[i].binary[j+n+1];
            z=z+guy[i-n-1].binary[j];
            h=h+guy[i+n+1].binary[j];
        }
    }

}*/
//there are 2 ways to find boundary, the second way is finding the most left-top point of the boundary, then move
//along the right side
void boundary(){
    //bull 用来跳出循环 xyzh用来判断边界
    int i,j,n,x,y,z,h,bull=0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            guy[i].boundary[j]=0;
            //printf("%d ",guy[i].boundary[j]);
        }
        //puts("");
    }
    //there are 2 ways to find boundary, the second way is find the most left-top point of the boundary, then
    //along the right side
    for (i=X; i<N-X; i++) {
        for (j=X; j<N-X; j++) {
            x=0;
            y=0;
            z=0;
            h=0;
            for (n=0; n<X; n++) {
                x=x+guy[i].binary[j-n-1];
                y=y+guy[i].binary[j+n+1];
                z=z+guy[i-n-1].binary[j];
                h=h+guy[i+n+1].binary[j];
            }
            if (x==0&&y==X&&guy[i].binary[j]==1) {
                guy[i].boundary[j]=1;
            }
            else if (y==0&&x==X&&guy[i].binary[j]==1){
                guy[i].boundary[j]=1;
            }
            else if (z==0&&h==X&&guy[i].binary[j]==1){
                guy[i].boundary[j]=1;
            }
            else if (h==0&&z==X&&guy[i].binary[j]==1){
                guy[i].boundary[j]=1;
            }
        }
    }
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            //printf("%d ",guy[i].boundary[j]);
            if (guy[i].boundary[j]==1) {
                bull=1;
                start.firsti=i;
                start.firstj=j;
                break;
            }
        }
        if (bull==1) {
            break;
        }
        //puts("");
    }
}
// the way to get outer boundary is that
// from left to right and right to left
// find out the first 1 and continue until 0
void outer_boundary(){
    int i,j,m[N][N]={0},bull=0,x=0,y=0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (guy[i].boundary[j]==1&&bull==0) {
                m[i][j]=guy[i].boundary[j];
                bull=1;
            }
            else if (guy[i].boundary[j]==1&&bull==1){
                m[i][j]=guy[i].boundary[j];
                bull=1;
            }
            else if (guy[i].boundary[j]==0&&bull==1){
                bull=0;
                break;
            }
        }

    }
    bull=0;
    for (i=0; i<N; i++) {
        for (j=N-1; j>=0; j--) {
            if (guy[i].boundary[j]==1&&bull==0) {
                m[i][j]=guy[i].boundary[j];
                bull=1;
            }
            else if (guy[i].boundary[j]==1&&bull==1){
                m[i][j]=guy[i].boundary[j];
                bull=1;
            }
            else if (guy[i].boundary[j]==0&&bull==1){
                bull=0;
                break;
            }
        }
        
    }
    
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            guy[i].boundary[j]=m[i][j];
            //printf("%d ",m[i][j]);
        }
        //puts("");
    }
    // add points at junction
    bull=0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            ///up right
            if (guy[i].boundary[j]==1&&(i!=x+1||j!=y-1)&&guy[i-1].boundary[j+1]==1&&bull==0) {
                guy[i].boundary[j+1]=1;
                x=i;
                y=j+1;
                bull=1;
            }
            //bot right
            else if (guy[i].boundary[j]==1&&guy[i+1].boundary[j+1]==1&&bull==0) {
                guy[i].boundary[j+1]=1;
                bull=1;
            }
            else if (guy[i].boundary[j]==0&&bull==1){
                bull=0;
                break;
            }
        }
    }
    bull=0;
    for (i=0; i<N; i++) {
        for (j=N-1; j>=0; j--) {
            ///up left
            if (guy[i].boundary[j]==1&&(i!=x+1||j!=y+1)&&guy[i-1].boundary[j-1]==1&&bull==0) {
                guy[i].boundary[j-1]=1;
                x=i;
                y=j-1;
                bull=1;
            }
            //bot left
            else if (guy[i].boundary[j]==1&&guy[i+1].boundary[j-1]==1&&bull==0) {
                guy[i].boundary[j-1]=1;
                bull=1;
            }
            else if (guy[i].boundary[j]==0&&bull==1){
                bull=0;
                break;
            }
        }
        
    }

}

void subsample(){
    // f 测试有多少点被取到
    int i,j,count=0,number,x;
    float d,standard = 570*570*2;
    number=N/B+1;
    number=pow(number,2);
    int temp[number][2];
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            guy[i].subsample[j]=0;
        }
    }
    //printf("%d\n",number);
    count=0;
    for (i=0; i<N; i=i+B) {
        for (j=0; j<N; j=j+B) {
            temp[count][0]=i;
            temp[count][1]=j;
            count++;
        }
    }
    //printf("%d\n",count);
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (guy[i].boundary[j]==1) {
                for (x=0; x<number; x++) {
                    d=fabsf(powf((i-temp[x][0]),2)+powf((j-temp[x][1]),2));
                    if (standard>d) {
                        standard=d;
                        start.selecti=temp[x][0];
                        start.selectj=temp[x][1];
                    }
                }
                standard = N*N*2;
            }
            if (start.selecti>0&&start.selectj>0) {
                guy[start.selecti].subsample[start.selectj]=1;
            }
            /*if ((start.selecti>0&&start.selectj>0)&&(start.selecti<550&&start.selectj<=550)) {
                guy[start.selecti].subsample[start.selectj]=1;
            }
            // the 500 350 point is coverd by some 550 xxx points
            // so I eliminate the 550 xxx points and add 500 350 point
            guy[500].subsample[350]=1;*/
        }
    }
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            //printf("%d ",guy[i].subsample[j]);
            if (guy[i].subsample[j]==1) {
                start.numofsub++;
                printf("%d %d\n",i,j);
            }
        }
        //puts("");
    }
    printf("%d %d %d",start.numofsub,start.firsti,start.firstj);
    puts("");
}
void writeExcel()
{
    int i,j ;
    FILE *fp = NULL ;
    fp = fopen(OUT,"w") ;
    for (i=0 ; i<N ;i++){
        for (j=0; j<N; j++) {
            fprintf(fp, "%d\t",guy[i].subsample[j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void chaincode(){
    int i,j,n=0,tempi=0,tempj=0,starti=0,startj=0,lasti=0,lastj=0;
    for (i=0; i<N; i=i+B) {
        for (j=0; j<N; j=j+B) {
            if (n==0&&guy[i].subsample[j]==1) {
                starti=i;
                startj=j;
                n=1;
            }
            if (n==1) {
                break;
            }
        }
        if (n==1) {
            n=0;
            break;
        }
    }
    n=0;
    tempi=starti;
    tempj=startj;
    while (1) {
        if (n!=0&&tempi==starti&&tempj==startj) {
            break;
        }
        if (guy[tempi].subsample[tempj]==1&&guy[tempi].subsample[tempj+50]==1&&(tempi!=lasti||tempj+50!=lastj)) {
            start.chaincode[n]=0;
            n++;
            lasti=tempi;
            lastj=tempj;
            tempi=tempi;
            tempj=tempj+50;
            //printf("%d ",start.chaincode[n-1]);
            continue;
        }
        else if(guy[tempi].subsample[tempj]==1&&guy[tempi+50].subsample[tempj]==1&&(tempi+50!=lasti||tempj!=lastj)){
            start.chaincode[n]=3;
            n++;
            lasti=tempi;
            lastj=tempj;
            tempi=tempi+50;
            tempj=tempj;
            //printf("%d ",start.chaincode[n-1]);
            continue;
        }
        else if(guy[tempi].subsample[tempj]==1&&guy[tempi].subsample[tempj-50]==1&&(tempi!=lasti||tempj-50!=lastj)){
            start.chaincode[n]=2;
            n++;
            lasti=tempi;
            lastj=tempj;
            tempi=tempi;
            tempj=tempj-50;
            //printf("%d ",start.chaincode[n-1]);
            continue;
        }
        else if(guy[tempi].subsample[tempj]==1&&guy[tempi-50].subsample[tempj]==1&&(tempi-50!=lasti||tempj!=lastj)){
            start.chaincode[n]=1;
            n++;
            lasti=tempi;
            lastj=tempj;
            tempi=tempi-50;
            tempj=tempj;
            //printf("%d ",start.chaincode[n-1]);
            continue;
        }
    }
    for (n=0; n<start.numofsub; n++) {
        printf("%d\t",start.chaincode[n]);
    }
}
