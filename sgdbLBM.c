#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//----------------------------------------------------------------------
#define nx 1100		  
#define ny 	 80			  
#define q		  9
//----------------------------------------------------------------------
const int time			= 200000;
const int noOfSnaps	= 		11;
const int dispFreq	= 	 200;
const double D			= 	20.0;  
const double tau		=		 0.8;  
const double rho0		= 	 1.0;  
const double u_w		=		0.05;  
const double xc			=	 200.0;  
const double yc			=		20.0;
//---------------------------------------------------------------------- 
int main(void)
{
	int i,j,a,a1,ts,ia,ja;  
	static int ex[q],ey[q],kb[q];  
	static int isn[nx+2][ny+2]; 

	double mass, MoI, tmp1, tmp2, tmp3, rhoAvg, fx_t, fy_t, torq_t, pi;  
	double x[2], y[2], u[3], v[3], omg[3], fx[2], fy[2], torq[2], ubx, uby;  
	static double f[q][nx+2][ny+2], ft[q][nx+2][ny+2], feq[q][nx+2][ny+2]; 
	static double wt[q], ux[nx+2][ny+2], uy[nx+2][ny+2], rho[nx+2][ny+2];  

	char prefix[]="snap_",type[]=".dat",filename[15],solstr[5];  
	int solnumber=0;  
	FILE *soln, *fid;
//----------------------------------------------------------------------
	ex[0] =	0; ey[0] = 0;
	ex[1] = 1; ey[1] = 0;
	ex[2] = 0; ey[2] = 1;
	ex[3] =-1; ey[3] = 0;
	ex[4] = 0; ey[4] =-1;
	ex[5] = 1; ey[5] = 1;
	ex[6] =-1; ey[6] = 1;
	ex[7] =-1; ey[7] =-1;
	ex[8] = 1; ey[8] =-1;
//----------------------------------------------------------------------
	for (a=0; a<9; a++) 
	{
		if (a==0) 				{wt[a] = 4.0/9.0 ;}
		if (a>=1 && a<=4)	{wt[a] = 1.0/9.0 ;}
		if (a>=5 && a<=8)	{wt[a] = 1.0/36.0;}
	}
//----------------------------------------------------------------------
  for (a=0; a<q; a++) 
	{
		for (a1=a; a1<q; a1++)
	{ 
			if ( ex[a]+ex[a1]==0 && ey[a]+ey[a1]==0)
			{
				kb[a]	= a1;
				kb[a1]=  a;
			}
		}
	}
//----------------------------------------------------------------------
	for (i=0; i<=nx+1; i++)
	{
		for (j=0; j<=ny+1; j++)
		{
			for (a=0; a<9; a++)
			{
				f[a][i][j] = wt[a]*rho0;
			}
		}
	}
//----------------------------------------------------------------------		
	x[0] = xc; 
	y[0] = yc;
	
	u[0] = 0.0; 
	v[0] = 0.0;
	u[1] = 0.0; 
	v[1] = 0.0;
	
	omg[0] = 0.0;
	omg[1] = 0.0;
	
	fx[0]  = 0.0;
	fy[0]  = 0.0;
	torq[0]= 0.0;	

	pi 		= acos(-1.0);
	mass 	= rho0*pi*D*D/4.0;
	MoI 	= mass*D*D/8.0;

	fid = fopen("tRhoxyuvfxfyOmgT.dat","w"); 
	fprintf(fid,"Variables=time,rho,x,y,u,v,fx,fy,omega,torque\n");
//----------------------------------------------------------------------
	for (ts=0; ts<=time; ts++)
	{	
		for (i=1; i<=nx; i++)
		{
			for (j=0; j<=ny+1; j++)
			{
				tmp1 = pow((pow(i-x[0],2.0) + pow(j-y[0],2.0)),0.5); 

				if (tmp1 <= 0.5*D)	isn[i][j] = 1;
				else if (j==0 )			isn[i][j] = 2; 
				else if (j==ny+1)		isn[i][j] = 3; 
				else								isn[i][j] = 0; 
			}
		}
//----------------------------------------------------------------------
		rhoAvg = 0.0;

		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				tmp1 = 0.0; 
				tmp2 = 0.0; 
				tmp3 = 0.0; 
				
				for (a=0; a<q; a++) 
				{
					tmp1 += f[a][i][j]; 
					tmp2 += f[a][i][j]*ex[a]; 
					tmp3 += f[a][i][j]*ey[a]; 
				}

				rho[i][j] = tmp1; 
				ux[i][j]	= tmp2/tmp1; 
				uy[i][j]	= tmp3/tmp1; 

				rhoAvg	 += tmp1;
			}
		}
		
		rhoAvg /=(nx*ny);
//----------------------------------------------------------------------
		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++) 
				{
					tmp1					= ux[i][j]*ex[a]  + uy[i][j]*ey[a] ; 
					tmp2					= ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j]; 
					feq[a][i][j]	= wt[a]*rho[i][j]*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2); //scalar can be used
					ft[a][i][j]		= f[a][i][j] - (f[a][i][j]-feq[a][i][j])/tau ; //collision
				}
			}
		}
//----------------------------------------------------------------------
		for (i=1; i<=nx; i++) //Streaming post-collision
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++)
				{
					ia = i+ex[a]; 
					ja = j+ey[a];
					
					if (ia<1 )	{ ia = nx;  }
					if (ia>nx)	{ ia = 1 ;  } 

					f[a][ia][ja] = ft[a][i][j];
				}
			}
		}
//----------------------------------------------------------------------
		fx[1]		=0.0;
		fy[1]		=0.0;
		torq[1]	=0.0;
 		
		for (i=1; i<=nx; i++) //BC
		{
			for (j=1; j<=ny; j++)
			{
				if (isn[i][j] == 0)
				{
					for (a=0; a<q; a++)
					{
						ia = i+ex[a]; 
						ja = j+ey[a]; 

						if (ia<1 )	{ia = nx;}
						if (ia>nx)	{ia = 1 ;}

						if (isn[ia][ja] == 2) //Bottom Wall
						{
							f[kb[a]][i][j] = f[a][ia][ja] + 6.0*wt[a]*0.5*u_w*ex[kb[a]]; 
						}

						if (isn[ia][ja] == 3) //Top Wall
						{
							f[kb[a]][i][j] = f[a][ia][ja] - 6.0*wt[a]*0.5*u_w*ex[kb[a]];
						}

						if (isn[ia][ja]==1) //Cylinder
						{
						  ubx = (u[1] - omg[1]*(j+0.5*ey[a]-y[1]));
							uby = (v[1] + omg[1]*(i+0.5*ex[a]-x[1]));
							
							f[kb[a]][i ][j ] = ft[a		 ][i ][j ] + 6.0*wt[a]*rho0*(ubx*ex[kb[a]] + uby*ey[kb[a]]);
							f[a    ][ia][ja] = ft[kb[a]][ia][ja] - 6.0*wt[a]*rho0*(ubx*ex[kb[a]] + uby*ey[kb[a]]);						
							
							tmp1 = ex[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j] - 6.0*wt[a]*rho0*(ubx*ex[a] + uby*ey[a]));						
							tmp2 = ey[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j] - 6.0*wt[a]*rho0*(ubx*ex[a] + uby*ey[a]));
							
							fx[1] += tmp1;
							fy[1] += tmp2;
							
							torq[1] += (i+0.5*ex[a]-x[1])*tmp2 - (j+0.5*ey[a]-y[1])*tmp1;
						}					
					}
				}
			}
		}
//----------------------------------------------------------------------
		fx_t 	= 0.5*(fx[0]  +  fx[1]);
		fy_t 	= 0.5*(fy[0]  +  fy[1]);
		torq_t= 0.5*(torq[0]+torq[1]);
		
		u[2]=u[0] + fx_t/mass; 
		v[2]=v[0] + fy_t/mass;
	
		u[1]= 0.5*(u[0] + u[2]);
		v[1]= 0.5*(v[0] + v[2]);
			
		x[1]= x[0] + u[2];
		if(x[1] + 0.5*D >= nx) x[1] = 0.5*D;
		
		y[1]= y[0] + v[2];
	
		omg[2] = omg[0] + torq_t/MoI;
		omg[1] = 0.5*(omg[0] + omg[2]);
//----------------------------------------------------------------------		
		if (ts % dispFreq == 0)
		{	
			fprintf(fid,"%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",\
							ts, rhoAvg, x[1], y[1]/ny, u[1]/u_w, v[1]/u_w, fx_t, fy_t, omg[1], torq_t);
			printf("%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",\
							ts, rhoAvg, x[1], y[1]/ny, u[1]/u_w, v[1]/u_w);
		}
//----------------------------------------------------------------------
		fx[0] = fx[1];
		fy[0] = fy[1];
		
		torq[0] = torq[1];
		
		x[0] = x[1];
		y[0] = y[1];
		
		u[0] = u[2];
		v[0] = v[2];
		
		omg[0] = omg[2];
//----------------------------------------------------------------------
		if(ts<=time && ts % (time/(noOfSnaps-1)) == 0)
		{
			solnumber++; 

			strcpy(filename,prefix);
			sprintf(solstr,"%d",solnumber); 
			strcat(filename,solstr); 
			strcat(filename,type); 
			soln = fopen(filename,"w"); 

			fprintf(soln,"Variables=x,y,u,v,rho,region\n"); 
			fprintf(soln,"Zone I= %d,J= %d\n\n",nx,ny); 

			for (j=1; j<=ny; j++)
			{
				for (i=1; i<=nx; i++)
				{
					fprintf(soln,"%d %d %12.8f %12.8f %12.8f %d\n",i,j,ux[i][j],uy[i][j],rho[i][j],isn[i][j]); 
				}
				fprintf(soln,"\n"); 
			}
			fclose(soln);
			printf("snap %d recorded at time %d\n",solnumber,ts);
		}
//----------------------------------------------------------------------
	}//Time loop Ends

	fclose(fid); 
	return 0; 
}
