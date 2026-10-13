/*************************************************************************
 file name: matutil.c
*************************************************************************/

#include "matutil.h"

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
  
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      b[i][j] = a[i][j];
    } /* for */
  } /* for */

  return;

} /* matCopy */

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
  
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
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
  
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
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
  
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < p; j++)
    {
      c[i][j] = 0.0;
      
      for (k = 0; k < n; k++)
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

      n - The number of rows and the number of columns in matrix a.
      
      a - A pointer to a pointer of an array of values.

      aInv - A pointer to a pointer to storage for the inverse of a.

      buffer - A pointer to a buffer to be used for temporary storage.
      
  Outputs:
  
      status - The returned status.
        0 - Failure due to a singular matrix.
        1 - Success.
  
*************************************************************************/
int matInverse(int n,double **a,double **aInv,double *buffer)
{
  int status;
  int i, j, position;
  int lda;  
  int ipriv[100]; 
  double work[100];
  int lwork = 100;
    
  // This is a square matrix.
  lda = n; 
    
  // Copy the matrix into the linear buffer.
  position  = 0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      buffer[position] = a[0][i + (n * j)];
      position++;
    } // for
  } // for

  // Perform factorization.
  dgetrf_(&n,&n,buffer,&lda,&ipriv,&status);

  if (status == 0)
  {
    // Invert.
    dgetri_(&n,buffer,&lda,&ipriv,&work,&lwork,&status);

    if (status == 0)
    {
      // Copy the linear buffer into the output matrix.
      position  = 0;
      for (i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
        {
          aInv[0][i + (n * j)] = buffer[position];
          position++;
        } // for
      } // for
    } // if
  } // if

  if (status == 0)
  {
    status = 1;
  } // if
  else
  {
    status = 0;
  } // else

  return (status);

} /* matInverse */
