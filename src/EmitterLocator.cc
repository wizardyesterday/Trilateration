#include <stdio.h>
#include <math.h>
#include "matutil.h"
#include "EmitterLocator.h"

/*************************************************************************

  Name: EmitterLocator
  
  Purpose: The purpose of this function is to servd as a constructer
  of an EmmiterLocator instance.
  
  Calling Sequence: EmitterLocator()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
EmitterLocator::EmitterLocator(int m, int spaceIsInFeet)
{
  int i;

  for (i = 0; i < 3; i++)
  {
    deltaXy[i] = deltaXym[i];
  } /* for */
  
  for (i = 0; i < MAX_RECEIVERS; i++)
  {
    deltaR[i] = deltaRm[i];
    c[i] = cm[i];
  } /* for */

  for (i = 0; i < 3; i++)
  {
    cT[i] = cTm[i];
    cTxc[i] = cTxcm[i];
    cTxcInv[i] = cTxcInvm[i];
    cTxdeltaR[i] = cTxdeltaRm[i];    
  } /* for */

  /* allocate storage for the vectors */
  r0i = new double[m];
  rn0i = new double[m];
  
  return;
  
} // EmitterLocator

/*************************************************************************

  Name: ~EmitterLocator
  
  Purpose: The purpose of this function is to serve as a destructer
  of an EmitterLocator instance.
  
  Calling Sequence: ~EmitterLocator()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
EmitterLocator::~EmitterLocator(void)
{

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  // free storage for the vectors.
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  if (r0i != NULL)
  {
    delete [] r0i;
    r0i = NULL;
  } // if

  if (rn0i != NULL)
  {
    delete [] rn0i;
    rn0i = NULL;
  } // if
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
 
} // ~EmitterLocator

/*************************************************************************

  Name: computePosition
  
  Purpose: The purpose of this function is to compute the position of
  an emitter.  A least squares method is used so that RF multipath,
  timing jitter, and receiver position errors can be dealt with.
    
  Calling Sequence: status = computePosition(m,xn,yn,v,x,y)
  
  Inputs:

      m - The number of time differences.
      
      xn - The initial x position estimate of the emitter.
      
      yn - The nominal y position best estimate of the emitter.
      
      v - An array of 3 vectors that represent the positions of receivers
      in two space and the time difference of arrival between a receiver
      and the reference receiver.  Note that v[0] represents the
      coordinates of the reference receiver.
      
      x - A pointer to storage of the returned x coordinate of the
      emitter.
      
      y - A pointer to storage of the returned y coordinate of the
      emitter.
                 
  Outputs:
  
      status - The success of the computation.
        0 - Success.
        -1 - Failure due to a singular max.
        -2 - Failure do to too many iterations.
  
*************************************************************************/
int EmitterLocator::computePosition(int m,double xn,double yn,struct vector3 *v,
                       double *x,double *y)
{
  int i, done, loopCount;
  double dx, dy, lamda;
  int status;

  /* used for while loop */
  done = 0;
  loopCount = 0;
  
  /* compute r0i */
  for (i = 0; i < m; i++)
  {
    if (spaceIsInFeet)
    {
      /* convert to ft */
      r0i[i] = (100000 * v[i].dt) / 101670;
    } /* if */
    else
    {
      /* convert to m */
      r0i[i] = (100000 * v[i].dt) / 333564;
    } /* else */
  } /* for */
  
  while (!done)
  {
    for (i = 0; i < m; i++)
    {
      /* compute differential direction cosine max */
      computeDirection(xn,v[i].x,v[0].x,yn,
                          v[i].y,v[0].y,&c[i][0],&c[i][1]);    
 
      /* compute nominal differential pseudo range */
      rn0i[i] = computeNominalRange(xn,v[i].x,v[0].x,yn,v[i].y,v[0].y); 

      /* compute incremental differential pseudo range */
      deltaR[i][0] = r0i[i] - rn0i[i];
    } /* for */
   
    /* solve the least squares equation for dx and dy */
    status = computeDeltaXy(m,c,deltaR,deltaXy);  

    /* a singular max will cause a status of zero */
    if (status == 0)
    { /* set value and bail out */
      *x = 0;
      *y = 0;
      /* indicate singular max */
      status = -1;
      /* bail out */
      break;
    } /* if */
    
    /* compute new nominal position location */
    xn += deltaXy[0][0];
    yn += deltaXy[1][0];

    /* compute distance increment */
    dx = deltaXy[0][0];
    dy = deltaXy[1][0];
    lamda = sqrt((dx * dx) + (dy * dy));    

    /* break out of loop when desired resolution is met */
    if (lamda < LAMDA_TOLERANCE)
    {
      /* set return values */
      *x = xn;
      *y = yn; 
      /* indicate success */
      status = 0;     
      /* bail out */
      break;
    } /* if */ 
     
    /* indicate that one more iteration occurred */
    loopCount++;
    
    /* break out of loop if too many iterations have occurred */
    if (loopCount >= MAX_ITERATIONS)
    {
      /* indicate too many iterations */
      status = -2;      
      /* bail out */
      break;
    } /* if */ 
  } /* while */
  
  return (status);
  
} /* computePosition */

/*************************************************************************

  Name: computeNomRange
  
  Purpose: The purpose of this function is to compute the nominal
  differential pseudo range.
  
  Calling Sequence: rn0i = computeNomRange(xn,xi,x0,yn,yi.y0)
  
  Inputs:

     xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.
           
  Outputs:
  
      rn0i - The nominal differential pseudo range.
  
*************************************************************************/
double EmitterLocator::computeNominalRange(double xn,double xi,double x0,
                             double yn,double yi,double y0)
{
  double xni, xn0;
  double yni, yn0;
  double dni, dn0;
  double rn0i;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0));
  
  /* compute pseudo range */
  rn0i = dni - dn0;
    
  return (rn0i);
  
} /* computeNomRange */

/*************************************************************************

  Name: computeDirection
  
  Purpose: The purpose of this function is to compute the differential
  cosine and differential sine associated with a particular receiver
  position with respect to the reference receiver and the emitter.
  
  Calling Sequence: computeDirection(xn,xi,x0,yn,yi,y0,cx,cy)
  
  Inputs:

      xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.

      cx - A pointer to the direction cosine associated with the receiver.
      
      cy - A pointer to the direction sine associated with the reciever.
  
                    
  Outputs:
  
      None.
  
*************************************************************************/
void EmitterLocator::computeDirection(double xn,double xi,double x0,
                         double yn,double yi,double y0,
			 double *cx,double *cy)
{
  double xni, xn0;
  double yni, yn0;
  double dni, dn0;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0));
  
  /* compute direction cosine */
  *cx = (xni / dni) - (xn0 / dn0);
  
  /* compute direction sine */
  *cy = (yni / dni) - (yn0 / dn0);
  
  return;
  
} /* computeDirection */

/*************************************************************************

  Name: computeDeltaXy
  
  Purpose: The purpose of this function is to compute the incremental
  position correction to a position estimate.
  
  Calling Sequence: status = computeDeltaXy(m,c,deltaR,DeltaXy)
  
  Inputs:

      m - The number of time differences.
      
      c - An array of pointers to the rows of c[1..m][1..2].
      
      deltaR - An array of pointers to the rows of deltaR[1..m][1..1].
      
      DeltaXy - An array of pointers to the rows of DeltaXy[1..2][1..1].
     
      
  Outputs:
  
      status - The success of the computation.
        0 - Failure due to a singular max.
        1 - Success.
  
*************************************************************************/
int EmitterLocator::computeDeltaXy(int m,double **c,double **deltaR,double **DeltaXy)
{
  int status;
  double buffer[4];
  double *aPtr, *bPtr;

  aPtr = &cTxcm[0][0];
  bPtr = &cTxcInvm[0][0];

  /* compute [C'] */
  matTranspose(m,2,c,cT);
  
  /* compute cTxc = [C'][C] */
  matMultiply(2,m,2,cT,c,cTxc);
  
  /* compute cInv = Inv([C'][C]) */
  status = matInverse(2,&aPtr,&bPtr,buffer);
  
  if (status == 1)
  {
    /* compute cTxdeltaR = [C'][deltaR] */
    matMultiply(2,m,1,cT,deltaR,cTxdeltaR);
  
    /* compute [DeltaXy] = Inv([C'][C])[C'][deltaR] */
    matMultiply(2,2,1,cTxcInv,cTxdeltaR,DeltaXy);
  } /* if */
  
  return (status);

} /* computeDeltaXy */
