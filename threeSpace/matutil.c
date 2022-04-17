/*************************************************************************
 file name: matutil.c
*************************************************************************/

#include "nr.h"
#include "nrutil.h"
#include "matutil.h"

static double amCopy[MAT_MAX_SIZE][MAT_MAX_SIZE];
static double *aCopy[MAT_MAX_SIZE]; 

/*************************************************************************

  Name: matInit
  
  Purpose: The purpose of this function is to initialize all static
  storage associated with the matrix system.
  
  Calling Sequence: matInit()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
void matInit(void)
{
  int i;

  for (i = 1; i < MAT_MAX_SIZE; i++)
  {
    aCopy[i] = amCopy[i];
  } /* for */
  
  return;
  
} /* matInit */

/*************************************************************************

  Name: matCopy
  
  Purpose: The purpose of this function is to copy the matrix a into b.
  
  Calling Sequence: matCopy(m,n,a,b)
  
  Inputs:

      m - The number of rows in matrix a and b.
      
      n - The number of columns in matrix a, and b.
      
      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..n], which
      represents the copy of a.
      
  Outputs:
  
      None.
  
*************************************************************************/
void matCopy(int m,int n,double **a,double **b)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      b[i][j] = a[i][j];
    } /* for */
  } /* for */

  return;

} /* matAdd */

/*************************************************************************

  Name: matTranspose
  
  Purpose: The purpose of this function is to compute the transpose of a
  matrix.  
  
  Calling Sequence: matTranspose(m,n,a,aT)
  
  Inputs:

      m - The number of rows in a.
      
      n - The number of columns in a.
      
      a - An array of pointers to the rows of a[1..m][1..n].
        
      aT - An array of pointers to the rows of aT[1..n][1..m], which
      represents the transpose of a.
       
  Outputs:

      None.
         
*************************************************************************/
void matTranspose(int m,int n,double **a,double **aT)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      aT[j][i] = a[i][j];
    } /* for */
  } /* for */
  
  return;

} /* matTranspose */

/*************************************************************************

  Name: matAdd
  
  Purpose: The purpose of this function is to add two matrices,
  a and b.
  
  Calling Sequence: matAdd(m,n,a,b,c)
  
  Inputs:

      m - The number of rows in matrix a and b.
      
      n - The number of columns in matrix a, and b.
      
      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..n].
      
      c - An array of pointers to the rows of c[1..n][1..n], which
      represents the result a + b.
       
  Outputs:
  
      None.
  
*************************************************************************/
void matAdd(int m,int n,double **a,double **b,double **c)
{
  int i, j;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= n; j++)
    {
      c[i][j] = a[i][j] + b[i][j];
    } /* for */
  } /* for */

  return;

} /* matAdd */

/*************************************************************************

  Name: matMultiply
  
  Purpose: The purpose of this function is to multiply two matrices,
  a and b.
  
  Calling Sequence: matMultiply(m,n,p,a,b,c)
  
  Inputs:

      m - The number of rows in matrix a.
      
      n - The number of columns in matrix a, and number of rows in
      matrix b.
      
      p - The number of columns in matrix b.

      a - An array of pointers to the rows of a[1..m][1..n].
      
      b - An array of pointers to the rows of b[1..n][1..p].
      
      c - An array of pointers to the rows of c[1..m][1..p], which
      represents the result a * b.
       
  Outputs:
  
      None.
  
*************************************************************************/
void matMultiply(int m,int n,int p,double **a,double **b,double **c)
{
  int i, j, k;
  
  for (i = 1; i <= m; i++)
  {
    for (j = 1; j <= p; j++)
    {
      c[i][j] = 0.0;
      
      for (k = 1; k <= n; k++)
      {
        c[i][j] += a[i][k] * b[k][j]; 
      } /* for */
    } /* for */ 
  } /* for */
  
  return;
  
} /* matMultiply */

/*************************************************************************

  Name: matInverse
  
  Purpose: The purpose of this function is to invert a matrix a.
  
  Calling Sequence: status = matInverse(n,a,aInv)
  
  Inputs:

      n - The number of rows and columns in matrix a.
      
      a - An array of pointers to the rows of a[1..n][1..n].
      
      aInv - An array of pointers to the rows of aInv[1..n][1..n], which
      represents the result of the inverse of a.
       
  Outputs:
  
      status - The returned status.
        0 - Failure due to a singular matrix.
        1 - Success.
  
*************************************************************************/
int matInverse(int n,double **a,double **aInv)
{
  double col[10], d;
  int i, j, indx[10];


  matCopy(n,n,a,aCopy);
  
  if (!ludcmp(aCopy,n,indx,&d))
  {
    return (0);
  } /* if */
  
  for(j = 1; j <= n; j++)
  {
    for (i = 1; i <= n; i++)
    {
      col[i] = 0.0;
    } /* for */
    
    col[j] = 1.0;
    lubksb(aCopy,n,indx,col);
    
    for (i = 1; i <= n; i++)
    {
      aInv[i][j] = col[i];
    } /* for */
  } /* for */
  
  return (1);
  
} /* matInverse */
