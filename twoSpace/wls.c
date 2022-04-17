#include <stdio.h>
#include <math.h>
#include "matutil.h"
#include "wls.h"

int spaceIsInFeet;

double dXym[3][2], dRm[WLS_MAX_RECEIVERS][2];
double *dXy[3], *dR[WLS_MAX_RECEIVERS];

double cm[WLS_MAX_RECEIVERS][3];
double *c[WLS_MAX_RECEIVERS];

double cTm[3][WLS_MAX_RECEIVERS], cTxcm[3][3], cTxcInvm[3][3], cTxdRm[3][2];
double *cT[3], *cTxc[3], *cTxcInv[3], *cTxdR[3];

/*************************************************************************

  Name: wls
  
  Purpose: The purpose of this program is to compute the position of
  an RF emitter, in 2 space, using the time difference of arrivial (TDOA)
  technique.  The position parameters have units of feet, and
  the time parameters have units of nanoseconds.
  
  Calling Sequence: wls [-f] < inputFile
  
  Inputs:
    
    An optional flag (the -f) which indicates that position parameters are
    in feet, otherwise, they are in units of meters.
       
    A text file with the following format:
      n
      x0 y0 t00
      x1 y1 t01
      x2 y2 t02
        .
	.
      xn-1 yn-1 t0n-1
      	     
      The tuples, [xi yi t0i]:i=0,1,2,...n-1 correspond to the locations
      of receivers, to the differences of time of arrival between the
      reference receiver, located at (x0,y0), and a receiver at location
      (xi,yi). Note that the receiver, located at (x0,y0) is designated
      as the reference receiver for which all time differences are computed.
      The tuple [xe ye] is the initial estimate of the position of the
      emitter.  The value n indicates the number of receivers that were
      used for the collection of measurements      
  
  Outputs:
  
     Printed output of the solution [x1 y1], which locates the RF emitter. 

*************************************************************************/
int main(int argc, char **argv)
{
  char buffer[80];
  struct vector3 v[WLS_MAX_RECEIVERS];  
  int i, n;
  double xe, ye;
  double x, y;
  double distance, theta;
  int status;
  
  if (argc == 2)
  {
    spaceIsInFeet = 1;
  } /* if */
  else
  {
    spaceIsInFeet = 0;
  } /* else */
  
  /* initialize the system */
  wlsInit();  

  /* read number of receivers */  
  fgets(buffer,sizeof(buffer),stdin);
  sscanf(buffer,"%d",&n);
  
  /* zero initial estimate */
  xe = 0;
  ye = 0;
      
  /* read receiver locations and time differences */
  for (i = 0; i < n; i++)
  {
    fgets(buffer,sizeof(buffer),stdin);
    sscanf(buffer,"%lf %lf %lf\n",&v[i].x,&v[i].y,&v[i].dt);

    /* accumulate for later averaging */
    xe += v[i].x;
    ye += v[i].y;
  } /* for */
      
  /* compute initial estimate */
  xe /= n;
  ye /= n;
      
  /* determine location of emitter */
  status = wlsComputePosition((n-1),xe,ye,v,&x,&y);
    
  switch (status)
  {
    case 0:
      /* output in cartesian coordinates */
      printf("\n[x y] = [%f %f]\n",x,y);

      /* compute distance */
      distance = sqrt((x*x) + (y*y));

      /* fix to avoid division by zero */
      if (x == 0)
      {
        x = 0.001;
      } /* if */
    
      /* compute bearing */
      theta = atan(y / x) * 180 / M_PI;
    
      /* output in polar coordinates */
      printf("[r theta] = [%f %f]\n",distance,theta);
      break;
	
    case -1:
      printf("singular matrix\n");
      break;
	
    case -2:
      printf("too many iterations\n");
        break;
	
    default:
      printf("status: %d\n",status);
      break;
  } /* case */

  exit(0);
  
} /* main */

/*************************************************************************

  Name: wlsInit
  
  Purpose: The purpose of this function is to initialize all data
  structures.
  
  Calling Sequence: wlsInit()
  
  Inputs:

      None.
      
  Outputs:
  
      None.
  
*************************************************************************/
void wlsInit(void)
{
  int i;

  /* initialize the matrix system */
  matInit();
  
  for (i = 0; i < 3; i++)
  {
    dXy[i] = dXym[i];
  } /* for */
  
  for (i = 0; i < WLS_MAX_RECEIVERS; i++)
  {
    dR[i] = dRm[i];
    c[i] = cm[i];
  } /* for */

  for (i = 0; i < 3; i++)
  {
    cT[i] = cTm[i];
    cTxc[i] = cTxcm[i];
    cTxcInv[i] = cTxcInvm[i];
    cTxdR[i] = cTxdRm[i];    
  } /* for */

  return;
  
} /* wlsInit */

/*************************************************************************

  Name: wlsComputePosition
  
  Purpose: The purpose of this function is to compute the position of
  an emitter.  A least squares method is used so that RF multipath,
  timing jitter, and receiver position errors can be dealt with.
    
  Calling Sequence: status = wlsComputePosition(m,xn,yn,v,x,y)
  
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
        -1 - Failure due to a singular matrix.
        -2 - Failure do to too many iterations.
  
*************************************************************************/
int wlsComputePosition(int m,double xn,double yn,struct vector3 *v,
                       double *x,double *y)
{
  int i, done, loopCount;
  double *r0i, *rn0i;
  double dx, dy, lamda;
  int status;

  /* allocate storage for the vectors */
  r0i = (double *)malloc((m + 1) * sizeof(double));
  rn0i = (double *)malloc((m + 1) * sizeof(double));
  
  /* used for while loop */
  done = 0;
  loopCount = 0;
  
  /* compute r0i */
  for (i = 1; i <= m; i++)
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
    for (i = 1; i <= m; i++)
    {
      /* compute differential direction cosine matrix */
      wlsComputeDirection(xn,v[i].x,v[0].x,yn,
                          v[i].y,v[0].y,&c[i][1],&c[i][2]);    
 
      /* compute nominal differential pseudo range */
      rn0i[i] = wlsComputeNomRange(xn,v[i].x,v[0].x,yn,v[i].y,v[0].y); 

      /* compute incremental differential pseudo range */
      dR[i][1] = r0i[i] - rn0i[i];
    } /* for */
   
    /* solve the least squares equation for dx and dy */
    status = wlsComputeDeltaXy(m,c,dR,dXy);  

    /* a singular matrix will cause a status of zero */
    if (status == 0)
    { /* set value and bail out */
      *x = 0;
      *y = 0;
      /* indicate singular matrix */
      status = -1;
      /* bail out */
      break;
    } /* if */
    
    /* compute new nominal position location */
    xn += dXy[1][1];
    yn += dXy[2][1];

    /* compute distance increment */
    dx = dXy[1][1];
    dy = dXy[2][1];
    lamda = sqrt((dx * dx) + (dy * dy));    

    /* break out of loop when desired resolution is met */
    if (lamda < WLS_LAMDA_TOLERANCE)
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
    if (loopCount >= WLS_MAX_ITERATIONS)
    {
      /* indicate too many iterations */
      status = -2;      
      /* bail out */
      break;
    } /* if */ 
  } /* while */
  
  /* free storage for the vectors */
  free(r0i);
  free(rn0i);
  
  return (status);
  
} /* wlsComputePosition */

/*************************************************************************

  Name: wlsComputeNomRange
  
  Purpose: The purpose of this function is to compute the nominal
  differential pseudo range.
  
  Calling Sequence: rn0i = wlsComputeNomRange(xn,xi,x0,yn,yi.y0)
  
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
double wlsComputeNomRange(double xn,double xi,double x0,
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
  
} /* wlsComputeNomRange */
				  
/*************************************************************************

  Name: wlsComputeDirection
  
  Purpose: The purpose of this function is to compute the differential
  cosine and differential sine associated with a particular receiver
  position with respect to the reference receiver and the emitter.
  
  Calling Sequence: wlsComputeDirection(xn,xi,x0,yn,yi,y0,cx,cy)
  
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
void wlsComputeDirection(double xn,double xi,double x0,
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
  
} /* wlsComputeDirection */

/*************************************************************************

  Name: wlsComputeDeltaXy
  
  Purpose: The purpose of this function is to compute the incremental
  position correction to a position estimate.
  
  Calling Sequence: status = wlsComputeDeltaXy(m,c,dR,dXy)
  
  Inputs:

      m - The number of time differences.
      
      c - An array of pointers to the rows of c[1..m][1..2].
      
      dR - An array of pointers to the rows of dR[1..m][1..1].
      
      dXy - An array of pointers to the rows of dXy[1..2][1..1].
     
      
  Outputs:
  
      status - The success of the computation.
        0 - Failure due to a singular matrix.
        1 - Success.
  
*************************************************************************/
int wlsComputeDeltaXy(int m,double **c,double **dR,double **dXy)
{
  int status;

  /* compute [C'] */
  matTranspose(m,2,c,cT);
  
  /* compute cTxc = [C'][C] */
  matMultiply(2,m,2,cT,c,cTxc);
  
  /* compute cInv = Inv([C'][C]) */
  status = matInverse(2,cTxc,cTxcInv);
  
  if (status == 1)
  {
    /* compute cTxdR = [C'][dR] */
    matMultiply(2,m,1,cT,dR,cTxdR);
  
    /* compute [dXy] = Inv([C'][C])[C'][dR] */
    matMultiply(2,2,1,cTxcInv,cTxdR,dXy);
  } /* if */
  
  return (status);

} /* wlsComputeDeltaXy */
