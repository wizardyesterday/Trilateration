#include <stdio.h>
#include <math.h>
#include "matutil.h"
#include "wls.h"

int spaceIsInFeet;

double dXyzm[4][2], dRm[WLS_MAX_RECEIVERS][2];
double *dXyz[4], *dR[WLS_MAX_RECEIVERS];

double cm[WLS_MAX_RECEIVERS][4];
double *c[WLS_MAX_RECEIVERS];

double cTm[4][WLS_MAX_RECEIVERS], cTxcm[4][4], cTxcInvm[4][4], cTxdRm[4][2];
double *cT[4], *cTxc[4], *cTxcInv[4], *cTxdR[4];

/*************************************************************************

  Name: wls3
  
  Purpose: The purpose of this program is to compute the position of
  an RF emitter, in 3 space, using the time difference of arrivial (TDOA)
  technique.  The position parameters have units of feet or meters, and
  the time parameters have units of nanoseconds.
  
  Calling Sequence: wls [-f] < inputFile
  
  Inputs:
    
    An optional flag (the -f) which indicates that position parameters are
    in feet, otherwise, they are in units of meters.
       
    A text file with the following format:
      n
      x0 y0 z0 t00
      x1 y1 z1 t01
      x2 y2 z2 t02
        .
	.
      xn-1 yn-1 zn-1 t0n-1
      	     
      The tuples, [xi yi zi t0i]:i=0,1,2,...n-1 correspond to the locations
      of receivers, to the differences of time of arrival between the
      reference receiver, located at (x0,y0,z0), and a receiver at location
      (xi,yi,zi). Note that the receiver, located at (x0,y0,z0) is designated
      as the reference receiver for which all time differences are computed.
      The value n indicates the number of receivers that were used for the
      collection of measurements      
  
  Outputs:
  
     Printed output of the solution [x y z], which locates the RF emitter. 

*************************************************************************/
int main(int argc, char **argv)
{
  char buffer[80];
  struct vector4 v[WLS_MAX_RECEIVERS];  
  int i, n;
  double xe, ye, ze;
  double x, y, z;
  double distance, phi, theta;
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
  wlsInit3();  

  /* read number of receivers */  
  fgets(buffer,sizeof(buffer),stdin);
  sscanf(buffer,"%d",&n);
  
  /* zero initial estimate */
  xe = 0;
  ye = 0;
  ze = 0;
    
  /* read receiver locations and time differences */
  for (i = 0; i < n; i++)
  {
    fgets(buffer,sizeof(buffer),stdin);
    sscanf(buffer,"%lf %lf %lf %lf\n",&v[i].x,&v[i].y,&v[i].z,&v[i].dt);
    
    /* accumulate for later averaging */
    xe += v[i].x;
    ye += v[i].y;
    ze += v[i].z;
  } /* for */
  
  /* compute initial estimate */
  xe /= n;
  ye /= n;
  ze /= n;  
      
  /* determine location of emitter */
  status = wlsComputePosition3((n-1),xe,ye,ze,v,&x,&y,&z);

  switch (status)
  {
    case 0:
      /* output in cartesian coordinates */
      printf("\n[x y z] = [%f %f %f]\n",x,y,z);

      /* compute distance */
      distance = sqrt((x*x) + (y*y) + (z*z));

      /* fix to avoid division by zero */
      if (x == 0)
      {
        x = 0.001;
      } /* if */
    
      /* compute bearing */
      phi = atan(y / x) * 180 / M_PI;
      theta = atan(z / x) * 180 / M_PI;
    
      /* output in polar coordinates */
      printf("[r phi theta] = [%f %f %f]\n",distance,phi,theta);
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
void wlsInit3(void)
{
  int i;

  /* initialize the matrix system */
  matInit();
  
  for (i = 0; i < 4; i++)
  {
    dXyz[i] = dXyzm[i];
  } /* for */
  
  for (i = 0; i < WLS_MAX_RECEIVERS; i++)
  {
    dR[i] = dRm[i];
    c[i] = cm[i];
  } /* for */

  for (i = 0; i < 4; i++)
  {
    cT[i] = cTm[i];
    cTxc[i] = cTxcm[i];
    cTxcInv[i] = cTxcInvm[i];
    cTxdR[i] = cTxdRm[i];    
  } /* for */

  return;
  
} /* wlsInit */

/*************************************************************************

  Name: wlsComputePosition3
  
  Purpose: The purpose of this function is to compute the position of
  an emitter.  A least squares method is used so that RF multipath,
  timing jitter, and receiver position errors can be dealt with.
    
  Calling Sequence: status = wlsComputePosition3(m,xn,yn,zn,v,x,y)
  
  Inputs:

      m - The number of time differences.
      
      xn - The initial x position estimate of the emitter.
      
      yn - The nominal y position best estimate of the emitter.
      
      zn - The nominal z position best estimate of the emitter.
      
      v - An array of 4 vectors that represent the positions of receivers
      in three space and the time difference of arrival between a receiver
      and the reference receiver.  Note that v[0] represents the
      coordinates of the reference receiver.
      
      x - A pointer to storage of the returned x coordinate of the
      emitter.
      
      y - A pointer to storage of the returned y coordinate of the
      emitter.
      
      z - A pointer to storage of the returned z coordinate of the
      emitter.
                 
  Outputs:
  
      status - The success of the computation.
        0 - Success.
        -1 - Failure due to a singular matrix.
        -2 - Failure do to too many iterations.
  
*************************************************************************/
int wlsComputePosition3(int m,double xn,double yn,double zn,
                        struct vector4 *v,
		        double *x,double *y,double *z)
{
  int i, done, loopCount;
  double *r0i, *rn0i;
  double dx, dy, dz, lamda;
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
      wlsComputeDirection3(xn,v[i].x,v[0].x,
                           yn,v[i].y,v[0].y,
                           zn,v[i].z,v[0].z,
			   &c[i][1],&c[i][2],&c[i][3]); 

      /* compute nominal differential pseudo range */
      rn0i[i] = wlsComputeNomRange3(xn,v[i].x,v[0].x,
                                    yn,v[i].y,v[0].y,
				    zn,v[i].z,v[0].z); 

      /* compute incremental differential pseudo range */
      dR[i][1] = r0i[i] - rn0i[i];
    } /* for */

    /* solve the least squares equation for dx and dy */
    status = wlsComputeDeltaXyz(m,c,dR,dXyz); 

    /* a singular matrix will cause a status of zero */
    if (status == 0)
    { /* set value and bail out */
      *x = 0;
      *y = 0;
      *z = 0;
      /* indicate singular matrix */
      status = -1;
      /* bail out */
      break;
    } /* if */
            
    /* compute new nominal position location */
    xn += dXyz[1][1];
    yn += dXyz[2][1];
    zn += dXyz[3][1];

    /* compute distance increment */
    dx = dXyz[1][1];
    dy = dXyz[2][1];
    dz = dXyz[3][1];
    lamda = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    
    /* break out of loop when desired resolution is met */
    if (lamda < WLS_LAMDA_TOLERANCE)
    {
      /* set return values */
      *x = xn;
      *y = yn;
      *z = zn;      
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
  
} /* wlsComputePosition3 */

/*************************************************************************

  Name: wlsComputeNomRange3
  
  Purpose: The purpose of this function is to compute the nominal
  differential pseudo range.
  
  Calling Sequence: rn0i = wlsComputeNomRange3(xn,xi,x0,yn,yi.y0,zn,zi,z0)
  
  Inputs:

     xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.
      
      zn - The nominal z position best estimate of the emitter.
      
      zi - The z position of receiver i.
      
      z0 - The z position of the reference receiver.
                
  Outputs:
  
      rn0i - The nominal differential pseudo range.
  
*************************************************************************/
double wlsComputeNomRange3(double xn,double xi,double x0,
                           double yn,double yi,double y0,
                           double zn,double zi,double z0)
{
  double xni, xn0;
  double yni, yn0;
  double zni, zn0;
  double dni, dn0;
  double rn0i;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  zni = zn - zi;
  zn0 = zn - z0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni) + (zni * zni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0) + (zn0 * zn0));
  
  /* compute pseudo range */
  rn0i = dni - dn0;
    
  return (rn0i);
  
} /* wlsComputeNomRange3 */
				  
/*************************************************************************

  Name: wlsComputeDirection3
  
  Purpose: The purpose of this function is to compute the differential
  cosine and differential sine associated with a particular receiver
  position with respect to the reference receiver and the emitter.
  
  Calling Sequence: wlsComputeDirection3(xn,xi,x0,yn,yi,y0,zn,zi,z0,cx,cy,cz)
  
  Inputs:

      xn - The nominal x position best estimate of the emitter.
      
      xi - The x position of receiver i.
      
      x0 - The x position of the reference receiver.
      
      yn - The nominal y position best estimate of the emitter.
      
      yi - The y position of receiver i.
      
      y0 - The y position of the reference receiver.

      zn - The nominal z position best estimate of the emitter.
      
      zi - The z position of receiver i.
      
      z0 - The z position of the reference receiver.

      cx - A pointer to the direction cosine/sine associated with the receiver.
      
      cy - A pointer to the direction cosine/sine associated with the receiver.
  
      cz - A pointer to the direction cosine/sine associated with the receiver.
  
                    
  Outputs:
  
      None.
  
*************************************************************************/
void wlsComputeDirection3(double xn,double xi,double x0,
                          double yn,double yi,double y0,
                          double zn,double zi,double z0,
			  double *cx,double *cy, double *cz)
{
  double xni, xn0;
  double yni, yn0;
  double zni, zn0;
  double dni, dn0;

  /* compute distance differences */  
  xni = xn - xi;
  xn0 = xn - x0;
  yni = yn - yi;
  yn0 = yn - y0;
  zni = zn - zi;
  zn0 = zn - z0;
  
  /* compute distance from receiver i */
  dni = sqrt((xni * xni) + (yni * yni) + (zni * zni));
  
  /* compute distance from reference receiver */
  dn0 = sqrt((xn0 * xn0) + (yn0 * yn0) + (zn0 * zn0));
  
  /* compute direction cosine/sine */
  *cx = (xni / dni) - (xn0 / dn0);
  
  /* compute direction cosine/sine */
  *cy = (yni / dni) - (yn0 / dn0);
  
  /* compute direction cosine/sine */
  *cz = (zni / dni) - (zn0 / dn0);
  
  return;
  
} /* wlsComputeDirection3 */

/*************************************************************************

  Name: wlsComputeDeltaXyz
  
  Purpose: The purpose of this function is to compute the incremental
  position correction to a position estimate.
  
  Calling Sequence: status = wlsComputeDeltaXyz(m,c,dR,dXy)
  
  Inputs:

      m - The number of time differences.
      
      c - An array of pointers to the rows of c[1..m][1..3].
      
      dR - An array of pointers to the rows of dR[1..m][1..1].
      
      dXyz - An array of pointers to the rows of dXyz[1..3][1..1].
     
      
  Outputs:
  
      status - The success of the computation.
        0 - Failure due to a singular matrix.
        1 - Success.
  
*************************************************************************/
int wlsComputeDeltaXyz(int m,double **c,double **dR,double **dXyz)
{
  int status;

  /* compute [C'] */
  matTranspose(m,3,c,cT);
  
  /* compute cTxc = [C'][C] */
  matMultiply(3,m,3,cT,c,cTxc);
  
  /* compute cInv = Inv([C'][C]) */
  status = matInverse(3,cTxc,cTxcInv);

  if (status == 1)
  {  
    /* compute cTxdR = [C'][dR] */
    matMultiply(3,m,1,cT,dR,cTxdR);
  
    /* compute [dXyz] = Inv([C'][C])[C'][dR] */
    matMultiply(3,3,1,cTxcInv,cTxdR,dXyz);
  } /* if */
  
  return (status);

} /* wlsComputeDeltaXyz */
