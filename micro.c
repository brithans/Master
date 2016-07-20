/******************************************************************
micro.c
Written by Kelsey Hoffman
program to calulate the amplification due to microlensing for WD and solar type

the equations are based on that in Sahu 2003, eq. 14

to compile
gcc -lm -o micro micro.c

date: June 27, 2012
last update: July 11, 2012
****************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"

/*********************** MAIN PROGRAM ******************************/
int main(int argc, char **argv)
{

  int i, // a counter
    j, // a counter
    nring, //number of rings used to divide star
    mintheta, 
    ntheta, // number of theta bins
    nstep; // number of orbital steps

  double
    rsource, // radius of source, solar units
    msource, // mass of source, solar units
    rlens, // radius of lens, solar units
    mlens, // mass of lens, solar units
    period, // period of system (days)
    p_sec, //period in seconds 
    b_min, //impact parameter
    orbsep, // orbital seperation, AU
    yoffset, // how close tp crossing the point source
    inclin, // inclination
    time, // time of orbit
    dtime, // time step
    phase,  // phase of orbit
    p_off, //phase offset
    tflux, // total unobscured flux
    flux, // with lensing
    test,
    x, // location of lens
    y, // location of lens
    z, // location of lens
    thres;  // threshold amplification

 char file_out[]="microtest.dat";

 FILE *output;

 /* void subroutines */ 
 double unflux();
 double lflux(); // routine to calculate the lensed flux

 // set defaults of system parameters, source solar type lens WD
 // thres = 1e-7; //threshold amplification, if less don't calculate
 p_off = pi/2.0; // phase offset, so transit doesn't start a zero time
 nstep = 8000;
 nring = 5000;
 mintheta = 2; //minimum number of theta bins

 rsource = 1.0;
 msource = 1.0;
 rlens = 0.01;
 mlens = 0.6;
 b_min = 0.5; 
 period = 9.125;

  // read in input values
 for(i = 1; i<argc; i++){
   if(argv[i][0]=='-'){
     switch(argv[i][1]){
     case 'b': // the impact parameter
       sscanf(argv[i+1], "%lf", &b_min);
       break;
     }
   }
 }
  /******** output files *********/
 // printf("The output file is: %s \n", file_out);
 output = fopen(file_out, "w");
 
 /* parameters of the system */
 p_sec = period*s_day;
 orbsep = pow((SQR(p_sec)*G*(msource+mlens)*m_sun)/(4.0*pi*pi), 1.0/3.0)/au;
 yoffset = b_min * rsource * r_sun;
 inclin = pi/2.0 - asin(yoffset/(orbsep*au));

 /**** get unlensed/untransited flux for normalization ****/
 tflux = unflux(rsource,nring,mintheta);
 //printf("Total flux: %g \n", tflux);

 /***** calculate the lensing contribution */ 
 dtime = p_sec/nstep;
 time = 0.0; // initialize 
 /* loop through orbital phase */ 
 for(i = 0; i<=nstep; i++){
   phase = 2*pi*(time/p_sec - (int) (time/p_sec))-p_off;
   if(phase > -0.25 & phase < 0.25){
     x = orbsep * sin(phase);
     y = orbsep * cos(phase) * cos(inclin);
     z = orbsep * cos(phase);
     /* only calculate for when the lens plane in infront of source, transit */
     if(z > 0.0){
       flux = lflux(x,y,z,mlens,rsource, rlens, nring, mintheta);
     }
     else (flux = tflux);
     //fprintf(output, "%lf %lf \n", phase, flux/tflux);
     printf("%.8lf %.10lf \n", phase, flux/tflux);
   }
   time = time + dtime;
 }
 
}

 /*************** Routine for Lensed Flux *************/
double lflux(double x, double y, double z, double mlens, double rsource, double rlens, int nring, int mintheta){
  int i, j, ntheta;
  double D, Rs, Re, l, alpha, r, drad, aring, dtheta, theta, 
    zarea, xp, yp, yl, ylplus, ylminus, aplus, aminus, amp, pflux, intens;

  /* for limb darkening */ 
  double nl1, nl2, nl3, nl4, u, limb;
  double yl2; // for testing

  nl1 = 0.5118;
  nl2 = 0.0525;
  nl3 = 0.4590;
  nl4 = -0.2727;

  Rs = sqrt(SQR(x) + SQR(y))*au;
  D = z*au;
  Re = sqrt(4.0*G*mlens*m_sun*D/SQR(c));
  /* check threshold ??????? *****************/
  l = y * au;
  alpha = asin(l/Rs);
  /* divide up source */
  pflux = 0.0;
  drad = rsource/(nring); //ring width
  for(i = 0; i < nring; i ++){
    r = (drad*(i+1)+drad*i)/2.0; // mid point of ring
    aring = pi*(SQR(drad*(i+1))-SQR(drad*i));
    ntheta = mintheta * (i+1);
    dtheta = 2.0*pi/(ntheta);
    for(j = 0; j < ntheta; j++){
      theta = (dtheta * (j+1) + dtheta*(j))/2.0;
      zarea = aring/ntheta;
      xp = r * cos (theta); // use xp and yp for limb darkening
      yp = r * sin (theta);
      
      /* limb darkening */ 
      limb = 0.0; // initialize
      u = sqrt(1 - (SQR(xp) + SQR(yp)));
      limb = 1 - (nl1*(1-pow(u,0.5)) + nl2*(1-u) + nl3*(1-pow(u,1.5))
		  + nl4*(1-SQR(u)));
  
      yl = sqrt(SQR(Rs) + SQR(r*r_sun) + 2*Rs*r*r_sun*cos(alpha + theta));
      // printf("yl : %g \n", yl);
      //yl = sqrt(SQR(x*r_au -xp*rsource*r_sun) + 
      //  SQR(y*r_au-yp*rsource*r_sun));
      // printf("Other yl-yl2: %g \n", yl-yl2);
      ylplus = 0.5*(yl + sqrt(SQR(yl) + 4*SQR(Re)));
      /* check if occulted */ 
      if(ylplus > rlens*r_sun)
	aplus =  (SQR(yl/Re) + 2.0)/(2.0*(yl/Re)*sqrt(SQR(yl/Re) + 4.0)) + 0.5;
      else aplus = 0.0;
      ylminus = 0.5*(yl - sqrt(SQR(yl) + 4*SQR(Re)));
      /* check if occulted */ 
      if(fabs(ylminus) > rlens*r_sun)
	aminus= (SQR(yl/Re) + 2.0)/(2.0*(yl/Re)*sqrt(SQR(yl/Re) + 4.0)) - 0.5;
      else aminus = 0.0;
      amp = aplus + aminus;
  
      intens = 1.0; // intensity of the element
	 
      pflux = pflux + intens*zarea*amp*limb;
      //	 pflux = pflux + intens*zarea*limb;
    }
  }
  return pflux;
}

/******* Routine for Unobscured Flux ***************/
double unflux(double rsource,int  nring,int mintheta){
  int i, j;
  double 
    flux,
    drad,
    r,
    theta,
    aring,
    ntheta,
    dtheta,
    zarea,
    xp,
    yp,
    intens;

  /* for limb darkening */ 
  double nl1, nl2, nl3, nl4, u, limb;


  // double G, c, m_sun, r_sun, pi, r_au;
  //G = 6.674e-11;
  //c = 2.99792458e10;
  //m_sun = 1.9891e33;
  //r_sun = 6.96265e10;
  //pi = acos(-1.0);
  //r_au = 1.4959787069e13;

  nl1 = 0.5118;
  nl2 = 0.0525;
  nl3 = 0.4590;
  nl4 = -0.2727;

     flux = 0.0; //initialize
     drad = rsource/(nring); //ring width
     for(i = 0; i < nring; i ++){
       r = (drad*(i+1)+drad*i)/2.0; // mid point of ring
       aring = pi*(SQR(drad*(i+1))-SQR(drad*i));
       ntheta = mintheta * (i+1);
       dtheta = 2.0*pi/ntheta;
       for(j = 0; j < ntheta; j++){
	 theta = (dtheta * (j+1) + dtheta*(j))/2.0;
	 zarea = aring/ntheta;
	 xp = r * cos (theta); // use xp and yp for limb darkening
	 yp = r * sin (theta); 
	 intens = 1.0; // intensity of the element
	 /* ADD IN LIMB DARKENING */

       /* limb darkening */ 
	 limb = 0.0; // initialize
	 u = sqrt(1 - (SQR(xp) + SQR(yp)));
	 limb = 1 - (nl1*(1-pow(u,0.5)) + nl2*(1-u) + nl3*(1-pow(u,1.5))
		     + nl4*(1-SQR(u)));

	 flux = flux + intens*zarea*limb;
       }
     }
     return flux;
}
