/**************NEW PROGRAM TO NUMERICALLY INTEGRATE FLUX***************/
//compile with gcc -lm -o micro micro.c
//store data in numintdata.dat

//Written by Britta Hansen

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include <stdlib.h>

float SG(float D_S,float D_LS,float M_b,float R_a, float R_b,float d, float ydist,float DENOM)
{
//printf("D_S = %f\tD_LS = %f\td = %f\tydist = %f\t DENOM = %f\n",D_S,D_LS,d,ydist,DENOM);
	float r;
	float theta;
	float xdist;//called l in paper

	xdist = sqrt(d*d - ydist*ydist);
	
	float D_L = D_S - D_LS;
	float alpha;
	float a = xdist/d;
	if((xdist!=0)&&(d!=0)){alpha = asin(a);}
	if((xdist==0)||(d==0)){alpha = 0;}
	float yp;
	float ym;

	
	float D = D_LS*D_L/D_S;
	float R_E = 2*sqrt(M_b*G*D)/c;
	
	float NUMERATORA =0;
	float NUMERATORB =0;
	float y;
	float Ap;
	float Am;
	float A;
	float dr = R_a/100;
printf("R_E = %f\tD = %f\n"R_e,D);
	if(d!=0)
	{
		for(r=0;r<R_a;r += dr)
		{
		Ap = 0;
		Am = 0;
		NUMERATORA = 0;
//printf("d = %f\tr = %f\td/r = %f\n",d,r,d/r);
			for (theta=0;theta<2*pi;theta +=0.05)
			{
        			y = sqrt(d*d + r*r + 2*d*r*cos(alpha + theta));
				yp = 0.5*(y+sqrt(y*y + 4*R_E*R_E));
				ym = 0.5*(y-sqrt(y*y + 4*R_E*R_E));

				if((yp>R_b)&&(y!=0))
				{
        				Ap = (pow(y/R_E,2) + 2)/((2*y/R_E)*sqrt(pow(y/R_E,2) + 4))+ 0.5;
				}
				if((sqrt(pow(ym,2))>R_b)&&(y!=0))
				{

					Am = (pow(y/R_E,2) + 2)/((2*y/R_E)*sqrt(pow(y/R_E,2) + 4)) - 0.5;
				}
				if(((yp<=R_b)&&(sqrt(pow(ym,2))<=R_b))||(y==0))
				{
					A=0;
				}
				A = Ap + Am;
				NUMERATORA += A*0.05;
			//	printf("NUMERATORA = %f\n",NUMERATORA);
			}	
		//printf("B=%f\n",NUMERATORB);
			NUMERATORB += NUMERATORA*dr;
		}
	}	

	if(d==0)
	{
		for(r=0;r<R_a;r += dr)
		{
			y=sqrt(r*r);
			Ap = (pow(y/R_E,2) + 2)/((2*y/R_E)*sqrt(pow(y/R_E,2) + 4))+ 0.5;
			Am = (pow(y/R_E,2) + 2)/((2*y/R_E)*sqrt(pow(y/R_E,2) + 4)) - 0.5;
			A = Ap + Am;
			NUMERATORA = A*2*pi;
			NUMERATORB += NUMERATORA*dr;
		}
	}

	A = NUMERATORB/DENOM;


	return(A);
}

int main(int argc, char *argv[])
{
	float A,Z,L_a,L_b;
	float M_a,M_b;//mass of the primary and secondary
	double phi,theta;//zenith and polar angles respectively
	double psi;//angle of rotation in it's own circular orbit.
	float mind,maxd;//minimum value of z == d =0
	float d;//pprojected distance between the two centres
	float R_a,R_b;//Primary and secondary radii
	float Flux; //Flux as a function of time
	float A_a,A_b;//Area of the first and second
	float v, w; //actual and observed rotational velocities
	float i;//inclination in radians
	float r;//radius of orbit
	int j,t,ttime;
	int Phase=1;//t will be in minutes
	int mint=0;
	float x,y,z;//new coordinates;
	float D_S;//distance from primary to us
	float micro;//value of microlensing added according to SG

	float DENOMINATORA =0;
	float DENOMINATORB =0;
	float DENOM=0;
/**************************
Default setting:*/
v = 0;
w = 0;
Flux = 0;
R_b = 0;
R_a = 0;
i = 0;
ttime = 10000;
psi=pi;
/***********************/

	for (j = 0; j<argc; j++)
	{
   		switch(argv[j][1])
		{
			case 'a': // radius of a
   			sscanf(argv[j+1],"%f",&R_a);
   			break;

   			case 'b': // radius of b called r_s in SG
    			sscanf(argv[j+1],"%f",&R_b);
    			break;

			case 'c'://mass of a
			sscanf(argv[j+1], "%f",&M_a);
			break;
	
			case 'd'://mass of b
			sscanf(argv[j+1],"%f",&M_b);
			break;

  			case 't': // total time
  			sscanf(argv[j+1], "%d", &ttime);
  			break;
				  
			case 'r': //radius of orbit called D_LS in SG
			sscanf(argv[j+1], "%f", &r);
			break;
				   
			case 'l'://luminosity of the primary
			sscanf(argv[j+1], "%f", &L_a);
			break;
				   
			case 'm'://luminosity of the secondary
			sscanf(argv[j+1], "%f", &L_b);
			break;

			case 'v'://rotational velocity
			sscanf(argv[j+1], "%f", &v);
			break;

			case 'w'://observed rotational velocity
			sscanf(argv[j+1], "%f", &w);
			break;

			case 'i'://inclination
			sscanf(argv[j+1], "%f", &i);
			break;

			case 'f'://distance from primary to us
			sscanf(argv[j+1], "%f", &D_S);
			break;
			
			case 'y'://min time
			sscanf(argv[j+1], "%d", &mint);
			break;
				   
			default:
			break;
   				
			case 'h': // HELP! 
			printf("Requires r,a,b,t,l,m,v,i\n");
			printf("Must be in solar radii, solar luminosities, solar masses, separations in AU, velocity in m/s\n");
  			fprintf(stderr, "Quick Helps: \n");
   			fprintf(stderr, "\n");
   			fprintf(stderr, " -r radius of orbit in AU\n");
   			fprintf(stderr, " -a radius_a in r_sun \n");
   			fprintf(stderr, " -b radius_b in r_sun \n");
			fprintf(stderr, " -c mass of a in m_sun\n");
			fprintf(stderr, " -d mass of c in m_sun\n");
   			fprintf(stderr, " -t total time in seconds \n");
			fprintf(stderr, " -l luminosity of the primary in l_sun\n");
			fprintf(stderr, " -m luminosity of the secondary in l_sun\n");
			fprintf(stderr, " -v actual rotational velocity in m/s\n");
			fprintf(stderr, " -w observed rotational velocity\n");
			fprintf(stderr, " -i inclination in radians\n");
			fprintf(stderr, " -f distance from source to us in AU\n");
   			fprintf(stderr, "\n");
   			break;
		}
	}

/***********************************************************************************
x-distance travelled in one second = d(1-cos(phi))
phi radians = (v_rot)*(1/d)(1 s)
w = v*sin(i);
if v==0 - then d = r
else if v!==0 then d = period*r*(1-cos(v/r))/4 //for circular orbit 

We're looking down x-axis, primary is at origin. Positive y to our right. Positive z is up. Psi is angle particle traces out in the circle in it's plane, with psi = 0 being the beginning of Phase 1.
***************************************************************************************/
//phi = pi/2 - i;
//theta = 3*pi/2;

/******************UNIT CONVERSION******************************************************/

R_a = R_a*r_sun;
R_b = R_b*r_sun;
M_a = M_a*m_sun;
M_b = M_b*m_sun;
r = au*r;
L_a = L_a*l_sun;
L_b = L_b*l_sun;
D_S = D_S*au;

DENOM = R_a*2*pi;

FILE *outfile;
char *mode= "w";
outfile = fopen("numintdata.dat",mode);
if (outfile ==NULL)
{
	fprintf(stderr,"Can't open output.txt\n");
	exit (1);
}
	A_a = pi*R_a*R_a;
	A_b = pi*R_b*R_b;
	const float firstz = z;//first value given of d.
	const float firstFlux = L_a/A_a + L_b/A_b;
	micro = 1;
	maxd = r - r/100;//+10
	mind = r*sin(i) + r/100;
psi = mint*v/r + pi;
printf("psi  = %f\n",psi);
if((psi>=pi)&&(psi<3*pi/2))
	{Phase = 1;}
if((psi>=pi*3/2)&&(psi<2*pi))
	{Phase = 2;}
if((psi>=pi*2)&&(psi<5*pi/2))
	{Phase = 3;}
if((psi>=5*pi/2)&&(psi<3*pi))
	{Phase = 4;}
printf("Intersection at %f seconds\n",pi*r/(v*2));
	for(t=mint;t<ttime;t++)
	{	
		psi = t*v/r + pi;
		x = r*cos(i)*sin(psi);
		y = r*cos(psi); //projected coordinates
		z = tan(i)*x;//plane of the orbit
		d = sqrt(z*z + y*y);//R_S
printf("x = %f\ty= %f\tz=%f\td = %f\n",x,y,z,d);
		if(Phase == 1)
		{
			if (R_a + R_b <= d)//separated
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro;
			fprintf(outfile,"%d\t%f\t%d\tseparated\n",t,Flux,Phase);
			}
			
			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro ;
			fprintf(outfile,"%d\t%lf\t%d\tintersecting\\n",t,Flux,Phase);
			}
			
			if ((d+R_b)<R_a)//inside
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro;
			fprintf(outfile,"%d\t%lf\t%d\tinside\tpsi = %f\t%f\t%f\n",t,Flux,Phase,psi,d,y);
			}
			if ((psi>=pi*3/2)&&(psi<2*pi)){
			Phase = 2;
			}
		}
		if(Phase == 2)
		{
			if ((d+R_b)<R_a)//inside
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro;
			fprintf(outfile,"%d\t%lf\t%d\tinside\t%f\t%f\n",t,Flux,Phase,d,y);
			}

			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro;
			fprintf(outfile,"%d\t%lf\t%d\tintersecting\n",t,Flux,Phase);
			}

			if (R_a + R_b <= d)//separated
			{
			micro = SG(D_S,r,M_b,R_a,R_b,d,sqrt(y*y),DENOM);
			Flux = micro;
			fprintf(outfile,"%d\t%f\t%d\tseparated\n",t,Flux,Phase);
			}	
			if ((y<-r*cos(2*v/r))||(y>r*cos(2*v/r))){//d>=maxd
			Phase = 3;
			}
		}
		if(Phase == 3)
		{
			if (R_a + R_b <= d)//separated
			{
			Flux = 1;
			fprintf(outfile,"%d\t%f\t%d\tseparated\tr=%f\td=%f\ty=%f\n",t,Flux,Phase,r,d,y);
			}

			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{micro = 1;
			A = AREA(R_a,R_b,d);
			Z = (A_b - A)/(A_b*A_b);
			Flux = ((L_a/A_a + L_b*Z)/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tintersecting\n",t,Flux,Phase);
			}
			
			if ((d+R_b)<R_a)//inside
			{
			Flux = (L_a/A_a/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tinside\n",t,Flux,Phase);
			}
			if ((y<r*(100000)*sin(v/(2*r)))&&(y>-r*(100000)*sin(v/(2*r)))){//d<=mind
			Phase = 4;
			}
		}
		if(Phase == 4)
		{
			if ((d+R_b)<R_a)//inside
			{
			Flux = (L_a/A_a/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tinside\n",t,Flux,Phase);
			}

			if ((R_a + R_b > d)&&((d+R_b)>=R_a))//intersecting
			{
			A = AREA(R_a,R_b,d);
			Z = (A_b - A)/(A_b*A_b);
			Flux = ((L_a/A_a + L_b*Z)/firstFlux);
			fprintf(outfile,"%d\t%lf\t%d\tintersecting\n",t,Flux,Phase);
			}

			if (R_a + R_b <= d)//separated
			{
			Flux = micro;
			fprintf(outfile,"%d\t%f\t%d\tseparated\n",t,Flux,Phase);
			}
			if ((y<-r*cos(2*v/r))||(y>r*cos(2*v/r))){
			Phase = 1;
			}
		}
	}
	fclose(outfile);
	
	return (0);
}
