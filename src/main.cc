#include <stdio.h>
#include <math.h>
#include "matutil.h"
#include "EmitterLocator.h"

int spaceIsInFeet;

double dXym[3][2], dRm[MAX_RECEIVERS][2];
double *dXy[3], *dR[MAX_RECEIVERS];

double cm[MAX_RECEIVERS][3];
double *c[MAX_RECEIVERS];

double cTm[3][MAX_RECEIVERS], cTxcm[2][2], cTxcInvm[2][2], cTxdRm[3][2];
double *cT[3], *cTxc[3], *cTxcInv[3], *cTxdR[3];

extern "C"
{
  extern void MAIN__(void);
};

void MAIN__(void) {}

/*************************************************************************

  Name: trilaterate
  
  Purpose: The purpose of this program is to compute the position of
  an RF emitter, in 2 space, using the time difference of arrivial (TDOA)
  technique.  The position parameters have units of feet, and
  the time parameters have units of nanoseconds.
  
  Calling Sequence: tri [-f] < inputFile
  
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
  struct vector3 v[MAX_RECEIVERS];  
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
  triInit();  

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
  status = computePosition((n-1),xe,ye,v,&x,&y);
    
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

  return (0);
  
} /* main */

