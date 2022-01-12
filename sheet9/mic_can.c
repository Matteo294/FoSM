#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

/* rescaling: all masses / m --> equals m = 1 in the relevant eqs.; all positions / sigma; all energies / epsilon; all velocities / sqrt(epsilon/m) 
 */

/* Initialize those factors as global variables
 */
const double m = 6.69e-26; // kg
const double epsilon = 1.65e-21; // J
const double sigma = 3.4e-10; // m
const double vscale_mic_can = 157.05;
const double k_B = 1.38e-23; // J/K


/* this is our data type for storing the information for the particles
 */
typedef struct
{
  double pos[3];
  double vel[3];
  double acc[3];
  double pot;
} 
particle;



/* auxiliary function to create a Gaussian random deviate
 */
double gaussian_rnd(void)
{
  double x, y;

  do
    {
      x = drand48();
      y = drand48();
    }
  while(x == 0);

    return sqrt(-2. * log(x)) * cos(2. * M_PI * y);

}


/* This function initializes our particle set.
 */
void initialize(particle * p, int nperdim, double boxsize, double temp, double epsilon_in_Kelvin)
{
  int i, j, k, n, l = 0;
  double spacing = boxsize/nperdim;
  double sigma_prime = sqrt(temp/epsilon_in_Kelvin); // sqrt(k_BT/m)/vscale = sqrt(k_BT * m / m *   epsilon)
  for(i = 0; i < nperdim; i++)
    for(j = 0; j < nperdim; j++)
      for(k = 0; k < nperdim; k++)
        {
	     
	      p[n].pos[0] = spacing * i; // naturally dimensionless bc boxsize is dimless in main func
	      p[n].pos[1] = spacing * j;
	      p[n].pos[2] = spacing * k;
	  
	      for (l = 0; l < 3; l++)
            {

              p[n].vel[l] = sigma_prime * gaussian_rnd();
            }
	  

          n++;
        }
}


/* This function updates the velocities by applying the accelerations for the given time interval.
 */
void kick(particle *p, int ntot, double dt)
{
  int i, k;

  for(i = 0; i < ntot; i++)
    for(k = 0; k < 3; k++)
        p[i].vel[k] += p[i].acc[k] * 0.5*dt;
}



/* This function drifts the particles with their velocities for the given time interval.
 * Afterwards, the particles are mapped periodically back to the box if needed.
 */
void drift(particle * p, int ntot, double boxsize, double dt)
{
  int i, k;

  for(i = 0; i < ntot; i++)
    for(k = 0; k < 3; k++)
      {

        p[i].pos[k] += p[i].vel[k] * dt;

        // map them escaped particles back into box; positions can go from 0 to boxsize
        if(p[i].pos[k] >= boxsize)
            p[i].pos[k] -= boxsize;
        if(p[i].pos[k] < 0)
            p[i].pos[k] += boxsize;

      }
    

}



/* This function calculates the potentials and forces for all particles. For simplicity,
 * we do this by going through all particle pairs i-j, and then adding the contributions both to i and j.
 */
void calc_forces(particle * p, int ntot, double boxsize, double rcut)
{
  int i, j, k;
  double rcut2 = rcut * rcut;
  double r2, r6, r12, dr[3], a, V;


  /* first, set all the accelerations and potentials to zero */
  for (i = 0; i < ntot; i++)
  {

    p[i].pot = 0;
    
    for (j = 0; j < 3; j++)
    {

      p[i].acc[j] = 0;

    }

  }

  /* sum over all distinct pairs */
  for(i = 0; i < ntot; i++)
    for(j = i + 1; j < ntot; j++)
      {
        for(k = 0, r2 = 0; k < 3; k++)
          {
       
            dr[k] = p[i].pos[k] - p[j].pos[k]; 
                 
            /* periodic mapping, i.e. we consider the box centered around the current particel
             * if the distance extends half the box size (in pos/neg direction) we substract/add
             * one boxsize so that the centered box picture is realized
             */
            if(dr[k] > 0.5 * boxsize)
                dr[k] -= boxsize;
            if(dr[k] < -0.5 * boxsize)
                dr[k] += boxsize;
              
                 r2 += pow(dr[k],2);
                  
          }

         r6 = pow(r2,3);
         r12 = pow(r6,2);
         if(r2 < rcut2) // potential stays zero if this cond isn't satisfied; quadratic formulation sufficient
          {

            V = 4 * (1/r12 - 1/r6);
            p[i].pot += V;
            p[j].pot += V;

          }
        	   
        for(k = 0; k < 3; k++)
          {
    
            a = 24 * dr[k] * (2/pow(r2,7) - 1/pow(r2,4)); // acc = -Nabla V
            p[i].acc[k] += a;
            p[j].acc[k] -= a; // a_ij = V(r_ij) per def, and a_ij = -a_ji

          }
      }
}


/* This function calculates the total kinetic and total potential energy, averaged per particle.
 * It also returns the instantanous kinetic temperature.
 */
void calc_energies(particle * p, int ntot, double epsilon_in_Kelvin, double *ekin, double *epot, double *temp)
{
  *ekin = 0; // need to set 0 bc otherwise takes old values in // why do we actually even take ekin etc as input? so that it can be pointed to also outside the function?
  *epot = 0;
  //*temp = 0; // newly ascribed, so safe without setting it zero before
  int i, j;
  for(i = 0; i < ntot; i++)
  {
    for(j = 0; j < 3; i++)
    {

      *ekin += 0.5 * pow(p[i].vel[j],2);      

    }

    *epot += 0.5 * p[i].pot; // 1/2 bc only once per pair

  }

  *ekin = *ekin/ntot;
  *epot = *epot/ntot; // normalized per particle
  *temp = 2./3. * *ekin * epsilon_in_Kelvin; // * epsilon bc Ekin is in dimless and for that temp is in kelvin // formula 9.13 in script 

}




/*
 * main driver routine
 */
int main()
{
  double epsilon_in_Kelvin = 120.0;     /* energy scale in Lennard Jones potential       */
  int N1d = 8;                          /* particles per dimension                       */
  int N = N1d * N1d * N1d;              /* total particle number                         */
  double mean_spacing = 5.0;            /* mean spacing in dimensionaless internal units */
  double boxsize = N1d * mean_spacing;  /* dimensionless boxsize                         */
  double target_temp_in_Kelvin = 80.0;  /* target temperature                            */
  int nsteps = 60000;                   /* number of steps to take                       */
  double dt = 0.01;                     /* timestep size                                 */

  double ekin, epot, temp;
  int step;

  /* allocate storage for our particles */
  particle *p = calloc(N , sizeof(particle));

  /* let's initialize the particles */
  initialize(p, N1d, boxsize, target_temp_in_Kelvin, epsilon_in_Kelvin);

  /* calculate the forces at t=0 */
  calc_forces(p, N, boxsize, 10.0);

  /* create an output file */
  char fname[20];
  sprintf(fname, "output_%d.txt", (int) target_temp_in_Kelvin);
  FILE *fd = fopen(fname, "w");

  /* measure energies at beginning, and output this to the file and screen */
  calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
  fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);
  printf("%6d   %10g   %10g   %10g       %10g\n", 0, ekin, epot, ekin + epot, temp);


  /* now we carry out time integration steps using the leapfrog */
  for(step = 0; step < nsteps; step++)
    {
      kick(p, N, dt);

      drift(p, N, boxsize, dt);

      calc_forces(p, N, boxsize, 10.0);

      kick(p, N, dt);


      /* measure energies and output this to the file and screen */
      calc_energies(p, N, epsilon_in_Kelvin, &ekin, &epot, &temp);
      fprintf(fd, "%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
      printf("%6d   %10g   %10g   %10g       %10g\n", step + 1, ekin, epot, ekin + epot, temp);
    }

  fclose(fd);

  return 0;
}
