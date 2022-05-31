//gcc -O3 main.c `pkg-config --libs gsl`

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#define REP 1000   //number of simulations
#define STEPS 200 //number of compute steps
#define H 0.01   //stepsize

//parameter for the damped dynamic mass - spring system
#define K 9000.0
#define M 450.0
#define D_MIN 800.0
#define D_MAX 1200.0


//DGL for the damped dynamic mass - spring system
int f(double t, const double y[], double dydt[], void *params)
{
	double k, d, m;
	
	k = ((double*)params)[0];
	d = ((double*)params)[1];
	m = ((double*)params)[2];

	dydt[0] = y[1];
	dydt[1] =  - d / m * y[1] - k / m * y[0];
	
	return GSL_SUCCESS;
}

//jacobian matrix for the damped dynamic mass - spring system
int jacobi(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
	double k, d, m;
	gsl_matrix *m_ptr;
	gsl_matrix_view dfdy_mat;
	
	k = ((double*)params)[0];
	d = ((double*)params)[1];
	m = ((double*)params)[2];
	
	dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);

    	m_ptr = &dfdy_mat.matrix;	/* m_ptr points to the matrix */

    	/* fill the Jacobian matrix as shown */
    	gsl_matrix_set (m_ptr, 0, 0, 0.0);	/* df[0]/dy[0] = 0 */
    	gsl_matrix_set (m_ptr, 0, 1, 1.0);	/* df[0]/dy[1] = 1 */
    	gsl_matrix_set (m_ptr, 1, 0, - k / m); /* df[1]/dy[0] */
    	gsl_matrix_set (m_ptr, 1, 1, - d / m);     /* df[1]/dy[1] */

    	/* set explicit t dependence of f[i] */
    	dfdt[0] = 0.0;
 	dfdt[1] = 0.0;


    	return GSL_SUCCESS;
}

//save function; save t and y in file
void save(char *file, double *t, double *y, int n)
{
	FILE *fd;
	int i;
	
	fd=fopen(file,"w");
	
	for(i=0; i<n; i++)
	{
		
		fprintf(fd,"%E ",t[i]);
		fprintf(fd,"%E ",y[i]);
		
		fprintf(fd,"\n");
	}
	fclose(fd);

}

int main()
{
	int dim = 2;
	int rep;
	int steps;
	double *tvec;
	double *y0vec;
	const gsl_odeiv2_step_type *type_ptr;
	gsl_odeiv2_step *step_ptr;
	gsl_odeiv2_evolve *evolve_ptr;
	gsl_odeiv2_system my_system;
	double y[2];
	double t;
	double params[3];
	
	//vector for the time and the mean motion
	tvec = calloc(sizeof(double), STEPS);
	y0vec = calloc(sizeof(double), STEPS);

	
	type_ptr = gsl_odeiv2_step_rk4; //RK4 Solver
	step_ptr = gsl_odeiv2_step_alloc(type_ptr, dim); //stepping function
	evolve_ptr = gsl_odeiv2_evolve_alloc(dim); //evolution function
	

	//parameter vector for the DGL
	params[0] = K;
	params[2] = M; 
	
	//DGL-System (jacobian is not used)
	my_system.function = f;	/* the right-hand-side functions dy[i]/dt */
    	my_system.jacobian = NULL;//jacobi;	/* the Jacobian df[i]/dy[j] */
    	my_system.dimension = dim;	/* number of diffeq's */
    	my_system.params = params;	/* parameters to pass to rhs and jacobian */
	
	for(rep=0; rep<REP; rep++)
	{
		//set the state vector motion=0 and velocity=0.1 
		y[0] = 0.0;
		y[1] = 0.1;
		//set time
		t = 0.0;
	
		//set damping factor
		params[1] = (D_MAX - D_MIN) * ((double)rand() / (double)RAND_MAX) + D_MIN;
		
		//first values
		tvec[0] = t;
		y0vec[0] += y[0];
		
		//sim. loop
		for(steps=1; steps < STEPS; steps++)
		{
			//ODE-Solver
			gsl_odeiv2_evolve_apply_fixed_step (evolve_ptr, NULL, step_ptr,
		                            &my_system, &t, H, y);

			//set result for time t
		    	tvec[steps] = t;
			y0vec[steps] += y[0];
		}
	
	}
	
	//compute the mean motion
	for(steps=0; steps < STEPS; steps++)
		y0vec[steps] /= (double)REP;

	//save results
	save("daten.dat", tvec, y0vec, STEPS);
	
	//final steps
	gsl_odeiv2_evolve_free (evolve_ptr);
    	gsl_odeiv2_step_free (step_ptr);
    	
	return EXIT_SUCCESS;
}
