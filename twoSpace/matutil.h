/*************************************************************************
 file name: matutil.h
*************************************************************************/

#ifndef __MATUTIL_H__
#define __MATUTIL_H__

/*************************************************************************
 defines    
*************************************************************************/
#define MAT_MAX_SIZE (50 + 1)

/*************************************************************************
 structure definitions     
*************************************************************************/     

/*************************************************************************
 function prototypes     
*************************************************************************/
void matInit(void);
void matCopy(int m,int n,double **a,double **b);
void matTranspose(int m,int n,double **a,double **aT);
void vecDotProduct(int n,double *v1,double *v2,double *v1DotV2);
void matAdd(int m,int n,double **a,double **b,double **c);
void matMultiply(int m,int n,int p,double **a,double **b,double **c);
int matInverse(int n,double **a,double **aInv);

#endif /* __MATUTIL_H__ */
