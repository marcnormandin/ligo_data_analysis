//// Condition Number calculation for the detector network matrix.
//  Network matrix values depend on the sky locations. Well-coditioned and ill-conditioned
//  nature of the network matrix (value of the condtion number) gives and idea of invertability of 
//  the matrix according the sky location. 
//  Large condition number -> large numerical errors in estimations.
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

//// Shihan Weerathunga 
double conditionnum_cal(
						const char* detId[], 
						double declination_input, double rightascension_input,
						double polarization_angle_input)
{
	//char* detId[] = { "H1", "H2"};

	int N = 4; //length(detId);
	gsl_vector* U_vec_input = gsl_vector_alloc( N );
	gsl_vector* V_vec_input = gsl_vector_alloc( N );

	for (int i = 0; i < N; i++) {
		double u, v, F_Plus, F_Cross;
		antennapattern(declination_input,
	 				rightascension_input,
	                polarization_angle_input,
	                detId(id),
	                &u, &v, &F_Plus, &F_cross);
		gsl_vector_set(U_vec_input, i, u);
		gsl_vector_set(V_vec_input, i, v);
	}

	UdotU_input=  dot(U_vec_input,U_vec_input);
	UdotV_input=  dot(U_vec_input,V_vec_input);
	VdotV_input=  dot(V_vec_input,V_vec_input);
	 
	A_input = UdotU_input; 
	B_input = UdotV_input; 
	C_input = VdotV_input;

	////  Network vector non-zero elements
	M(1,1)= A_input ;         M(1,2)= B_input ;
	M(2,1)= B_input ;         M(2,2)= C_input ;
	M(3,3)= A_input ;         M(3,4)= B_input ;
	M(4,3)= B_input;          M(4,4)= C_input;

	cond_num = cond(M);
	return cond_num;
}
