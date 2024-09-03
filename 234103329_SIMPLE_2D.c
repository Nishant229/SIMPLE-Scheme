#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

double general_boundary_conditions( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2] );

double Initialization(int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N]);

double u_velocity_star_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N]);

double v_velocity_star_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N] );

double pressure_dash_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_dash[M][N] );

double Results_Collocated_Grid(int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity[M+1][N+1][2]);

double Simple_Algorithm( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure[M][N] );

double max(double a, double b, double c);
// For Testing Purpose
void print_u_velocity(int M, int N, double vel[M+1][N+1][2]);
void print_v_velocity(int M, int N, double vel[M+1][N+1][2]);
void print_pressure(int M, int N, double pres[M][N]);



int main() {
	//// Dimensions
	float L=1.0, H=1.0; // length of square cavity
	
	//// Taking input from user for Grid size
	int M, N;
	/* printf("Please enter the Number of Grid Points on X-axis: ");
	scanf("%d",&M);
	printf("Please enter the Number of Grid Points on Y-axis: ");
	scanf("%d",&N);
	*/
	M=128, N=128;
	
	//// Grid Size
	double dx=L/(M), dy=H/(N);
	double grid_size[2];
	grid_size[0] = dx;
	grid_size[1] = dy;

	printf("dx = %lf ,  dy = %lf  \n\n", dx,dy);
	
	//// Boundary Conditions
	double U_top=1, U_bottom=0, U_right=0, U_left=0; // For u velocity
	double V_top=0, V_bottom=0, V_right=0, V_left=0; // For v velocity
	// Boundary Conditions Matrix
	double BCs[][4] = {{U_top, U_bottom, U_right, U_left}, // For u velocity
	{V_top, V_bottom, V_right, V_left}}; // For v velocity
	
	
	//// Reynolds Number
	double Re;
	/*printf("Please enter the Reynold's Number for the flow: ");
	scanf("%lf",&Re);*/
	Re = 1000.0;
	
	//// Pressure Matrix
	double pressure[M][N];
	//// Velocity Matrix
	double velocity[M+1][N+1][2]; 
	// U_Veloxity[M+1][N][0] ; V_Veloxity[M][N+1][1] ;
	double velocity_coefficients[M][N][2]; // 6 coefficients -> ((a_point, a_east, a_west, a_north, a_south, source)) for a velocity cell center in 2D space
	// Only for storing point velocity which is used in pressure correction equation
	
	// Implementing Simple Algorithm
	int iter =  Simple_Algorithm( M, N, Re, grid_size, BCs, velocity, velocity_coefficients, pressure );
	
	// Printing Results on Collocated Grid
	Results_Collocated_Grid(M, N, Re, grid_size, BCs, velocity);
	
	
	return 1;
}

double Simple_Algorithm( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure[M][N] ) {
	
	double velocity_old[M+1][N+1][2], velocity_new[M+1][N+1][2];
	double pressure_old[M][N], pressure_new[M][N];
	double velocity_star[M+1][N+1][2];
	double pressure_star[M][N];
	double pressure_dash[M][N];
	double epsilon;
	if (Re==100 || Re==400 ) {
		epsilon = pow(10,-5);
	}
	if (Re==1000 ) {
		epsilon = pow(10,-5);
	}
	double error = 1.0, e=1.0;
	double p_error = 0.0, u_error = 0.0, v_error = 0.0;  // L2_norm
	
	double dx= grid_size[0];
	double dy= grid_size[1];
	
	int SIMPLE_iteration_count=0;
	
	// Step 0
	Initialization( M, N, Re, grid_size, BCs, velocity_star, velocity_coefficients, pressure_star);
	//printf("-------------Initialization Complete-------------\n");
	
	general_boundary_conditions( M, N, Re, grid_size, BCs, velocity_star, velocity_coefficients);
	//printf("-------------Boundary Condition Applied-------------\n");
		
	// Testing
	//printf("The Parameters before Simple Algorithm is as follows: ");
	//print_u_velocity(M,N,velocity_star);
	//print_v_velocity(M,N,velocity_star);
	//print_pressure(M,N,pressure_star)
	 
	while /*(SIMPLE_iteration_count<1) {*/(error>epsilon) {
		
		SIMPLE_iteration_count++;
		//printf("SIMPLE Algorithm's Iteration - %d ; Error = %lf\n",SIMPLE_iteration_count,error);
		printf("SIMPLE Algorithm's Iteration - %d ; p_error = %lf ; u_error = %lf ; v_error = %lf ; Error = %lf\n",SIMPLE_iteration_count,p_error,u_error,v_error, error);
		p_error = 0, u_error = 0, v_error = 0;
		error = 0;
		
		// Step 1
		u_velocity_star_update( M, N, Re, grid_size, BCs, velocity_star, velocity_coefficients, pressure_star);
		v_velocity_star_update( M, N, Re, grid_size, BCs, velocity_star, velocity_coefficients, pressure_star);
		//printf("-------------Step 1 Complete-------------\n");
		
		// Testing
		//printf("The Velocity_star after solving for velocity* is as follows: ");
		//print_u_velocity(M,N,velocity_star);
		//print_v_velocity(M,N,velocity_star);
		
		// Step 2
		// Initializing p_dash = 0
		for(int i=0;i<M;i++) {
			for(int j=0;j<N;j++) {
					pressure_dash[i][j] = 0;
			}
		}
		p_error = pressure_dash_update( M, N, Re, grid_size, BCs, velocity_star, velocity_coefficients, pressure_dash);
		//printf("-------------Step 2 Complete-------------\n");
		// Testing
		//printf("The Pressure after solving pressure correction equation is as follows: ");
		//print_pressure(M,N,pressure_dash);
		
		// step 3
		// Correcting Pressure and Velocity
		double alpha_p=1.0, alpha_u=1.0, alpha_v=1.0;
		if (Re==100 || Re==400) {
			alpha_p=0.5, alpha_u=0.5, alpha_v=0.5;
		}
		else if(Re==1000) {
			alpha_p=0.05, alpha_u=0.05, alpha_v=0.05;
		}
		for(int i=0;i<M;i++) {
			for(int j=0;j<N;j++) {
					// pressure correction
					pressure_new[i][j] = pressure_star[i][j] + ( alpha_p * pressure_dash[i][j] );
					//p_error = p_error + pow(pressure_new[i][j] - pressure_old[i][j] , 2);
					// u_velocity correction
					if(i>0 && i<M) {
						double d_i_J= dy / velocity_coefficients[i][j][0];
						velocity_new[i][j][0] = velocity_star[i][j][0] + ( d_i_J*(pressure_dash[i-1][j] - pressure_dash[i][j]) );
						velocity_new[i][j][0] = alpha_u*velocity_new[i][j][0] + (1-alpha_u)*velocity_old[i][j][0];
						u_error = u_error + pow(velocity_new[i][j][0]- velocity_old[i][j][0] , 2);
					}
					// v_velocity correction
					if(j>0 && j<N) {
						double d_I_j= dx / velocity_coefficients[i][j][1];
						velocity_new[i][j][1] = velocity_star[i][j][1] + ( d_I_j*(pressure_dash[i][j-1] - pressure_dash[i][j]) );
						velocity_new[i][j][1] = alpha_v*velocity_new[i][j][1] + (1-alpha_v)*velocity_old[i][j][1];
						v_error = v_error + pow(velocity_new[i][j][1]- velocity_old[i][j][1] , 2);
					}
			}
		}
		// Error Calculation
		p_error = pow( ( p_error / (M*N) ), 0.5);
		u_error = pow( ( u_error / ((M-1)*N) ), 0.5);
		v_error = pow( ( v_error / (M*(N-1)) ), 0.5);
		error = p_error;
		//if (Re==100 || Re==400) {
			if (u_error>error) { error = u_error; }
			if (v_error>error) { error = v_error; }
		//}
		// printf("-------------Step 3 Complete-------%lf------\n",error);
		
		// step 4 - Solving discretized Equation
		// step 5 - Convergence check
		
		//print_u_velocity(M,N,velocity_new);
		//print_v_velocity(M,N,velocity_new);
		// Testing
		/*if (SIMPLE_iteration_count == 7975) {
			print_u_velocity(M,N,velocity_star);
			print_v_velocity(M,N,velocity_star);
			print_pressure(M,N,pressure_dash);
		}*/
		
		// step 6 - Setting parameters for next iteration
		for(int i=0;i<M;i++) {
			for(int j=0;j<N;j++) {
					// pressure
					pressure_old[M][N] = pressure_new[M][N];
					pressure_star[i][j] = pressure_new[i][j];
					// u_velocity
					if(i>0) {
						velocity_old[i][j][0] = velocity_new[i][j][0];
						velocity_star[i][j][0] = velocity_new[i][j][0];
					}
					// v_velocity
					if (j>0) {
						velocity_old[i][j][1] = velocity_new[i][j][1];
						velocity_star[i][j][1] = velocity_new[i][j][1];
					}
			}
		}
		// Testing
		//printf("--------After %d Simple Iteration --------u & v velocity-----\n", SIMPLE_iteration_count);
		//print_u_velocity(M,N,velocity_star);
		//print_v_velocity(M,N,velocity_star);
		
	}
	// Testing
	//print_u_velocity(M,N,velocity_star);
	//print_v_velocity(M,N,velocity_star);
	//print_pressure(M,N,pressure_dash);
	
	general_boundary_conditions( M, N, Re, grid_size, BCs, velocity, velocity_coefficients);
	// Storing the Results
	for(int i=0;i<=M;i++) {
		for(int j=0;j<=N;j++) {
				if (j<N && i>0 && i<M) {
					velocity[i][j][0] = velocity_new[i][j][0]; // final u velocity 
				}
				if (i<M && j>0 && j<N) {
					velocity[i][j][1] = velocity_new[i][j][1]; // final v velocity 
				}
				if (i<M && j<N) {
					pressure[i][j] = pressure_new[i][j]; // final pressure
				}
		}
	}
	
	return SIMPLE_iteration_count;
}


double general_boundary_conditions( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2] ) {
	// Left and Right Boundary
	for(int j=0;j<N;j++) {
			velocity_star[0][j][0] = BCs[0][3]; // Left Boundary
			velocity_star[M][j][0] = BCs[0][2]; // Right Boundary
	}
	// Top and Bottom Boundary
	for(int i=0;i<M;i++) {
			velocity_star[i][N][1] = BCs[1][0]; // Top Boundary
			velocity_star[i][0][1] = BCs[1][1]; // Bottom Boundary
	}
	return 1;
}


double Initialization( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N] ) {
	
	// Initializing velocity and pressure
	for(int i=0;i<=M;i++) {
		for(int j=0;j<=N;j++) {
				if (j<=N ) {
					velocity_star[i][j][0] = 0; // u velocity initialization
				}
				if (i<=M ) {
					velocity_star[i][j][1] = 0; // v velocity initialization
				}
				if (i<=M && j<=N ) {
					pressure_star[i][j] = 0; // pressure initialization
				}
		}
	}
	
	return 1;
}

double u_velocity_star_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N]  ) {
	double dx= grid_size[0];
	double dy= grid_size[1];
	double velocity_star_old[M+1][N+1][2];
	for(int i=0;i<=M;i++) {
		for(int j=0;j<=N;j++) {
			velocity_star_old[i][j][0] = velocity_star[i][j][0];
			velocity_star_old[i][j][1] = velocity_star[i][j][1]; 
		}
	}
	
	// Initializing u velocity equation coefficients
	double Fe=0, Fw=0, Fn=0, Fs=0; // Fe= east face, Fw= west face, Fn= north face, Fs= south face
	double De=0, Dw=0, Dn=0, Ds=0; // De= east face, Dw= west face, Dn= north face, Ds= south face
	double ap=0, ae=0, aw=0, an=0, as=0, source=0; // ap=point, ae= east face, aw= west face, an= north face, as= south face, source = source term
	
	//// ** Calculation only for internal nodes **
	//// U velocity --> Exclude all boundary velocities from all 2 dimensions
	
	double Su = 0, Sp = 0;
	double vel=0;
	double w=1.0; // under_relaxation factor
	
	// Testing
	//printf("u velocity from last iteration before updating coefficients \n");
	//print_u_velocity(M, N, velocity_star);
	
	// U velocity
	for(int i=1;i<M;i++) {
		for(int j=0;j<N;j++) {
			// neighbouring coefficients parameters
					// Convective parameters
					Su = 0;
					Sp = 0;
					
					Fe=( velocity_star_old[i+1][j][0] + velocity_star_old[i][j][0] )*dy / 2.0;
					Fw=( velocity_star_old[i][j][0] + velocity_star_old[i-1][j][0] )*dy / 2.0;
					if (j==N-1) { // Top Boundary
						Fn=( BCs[1][0] )*dx; 
					} 
					else { 
						Fn=( velocity_star_old[i][j+1][1] + velocity_star_old[i-1][j+1][1] )*dx / 2.0; 
					}
					if (j==0) { // Bottom Boundary
						Fs=( BCs[1][1] )*dx; 
					} 
					else { 
						Fs=( velocity_star_old[i][j][1] + velocity_star_old[i-1][j][1] )*dx / 2.0; 
					}
					// Diffusive parameters
					De= dy / ( Re * dx);
					Dw= dy / ( Re * dx);
					if (j==N-1) { // Top Boundary
						Sp = Sp - 2.0*dx / (Re * dy);
						Su = Su + (2.0*BCs[0][0]*dx) / (Re * dy);
					}
					else { 
						Dn= 4*dx / ( 4 * Re * dy); 
					}
					if (j==0) { // Bottom Boundary
						Sp = Sp - 2.0*dx / (Re * dy);
						Su = Su + (2.0*BCs[0][1]*dx) / (Re * dy);
					}
					else { 
						Ds= 4*dx / ( 4 * Re * dy); 
					}
					// neighbouring coefficients
					ae = max(-Fe,(De-Fe/2.0),0);
					aw = max(Fw,(Dw+Fw/2.0),0);
					if (j==N-1) { // Top Boundary
						an=0;
					}
					else {
						an = max(-Fn,(Dn-Fn/2.0),0);
					}
					if (j==0) { // Bottom Boundary
						as=0;
					}
					else {
						as = max(Fs,(Ds+Fs/2.0),0);
					}
					// Point term
					velocity_coefficients[i][j][0] = ae + aw + an + as + ( Fe - Fw + Fn - Fs  ) - Sp;
					ap = ae + aw + an + as + ( Fe - Fw + Fn - Fs  ) - Sp;
					// Pressure term coefficients or source term
					source = dy*(pressure_star[i-1][j] - pressure_star[i][j]) + Su;
					
					//vel = (ae*velocity_star[i+1][j][0] + aw*velocity_star[i-1][j][0]  +
					//      an*velocity_star[i][j+1][0] + as*velocity_star[i][j-1][0] +   source) /   ap;
					
					// Solved using Gauss Seidel method
					vel = source;
					if(ae!=0) {vel = vel + ae*velocity_star[i+1][j][0];} 
					if(aw!=0) {vel = vel + aw*velocity_star[i-1][j][0];}
					if(an!=0) {vel = vel + an*velocity_star[i][j+1][0];}
					if(as!=0) {vel = vel + as*velocity_star[i][j-1][0];}   
					vel = vel / ap;
					velocity_star[i][j][0] = (1-w)*velocity_star_old[i][j][0] + w*vel;
					
					// Testing
					//printf("u_velocity coefficients ");
					//printf("(%d,%d) --> %lf --> ",i,j,velocity_star[i][j][0]);
					//printf("Fe=%lf , Fw=%lf , Fn=%lf , Fs=%lf   \n",Fe, Fw,Fn,Fs);
					//printf("ap=%lf , ae=%lf , aw=%lf , an=%lf , as=%lf , source=%lf \n",velocity_coefficients[i][j][0], ae,aw,an,as,source);
					//printf("ap=%lf \n",velocity_coefficients[i][j][0][0]);
		}
	}
	return 1;
}

double v_velocity_star_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_star[M][N] ) {
	double dx= grid_size[0];
	double dy= grid_size[1];
	double velocity_star_old[M+1][N+1][2];
	for(int i=0;i<=M;i++) {
		for(int j=0;j<=N;j++) {
			velocity_star_old[i][j][0] = velocity_star[i][j][0];
			velocity_star_old[i][j][1] = velocity_star[i][j][1]; 
		}
	}
	
	// Initializing v velocity equation coefficients
	double Fe=0, Fw=0, Fn=0, Fs=0; // Fe= east face, Fw= west face, Fn= north face, Fs= south face
	double De=0, Dw=0, Dn=0, Ds=0; // De= east face, Dw= west face, Dn= north face, Ds= south face
	double ap=0, ae=0, aw=0, an=0, as=0, source=0; // ap=point value, ae= east face, aw= west face, an= north face, as= south face, source = source value
	
	//// ** Calculation only for internal nodes **
	//// V velocity --> Exclude all boundary velocities from all 2 dimensions
	
	double Su = 0, Sp = 0;
	double vel=0;
	double w=1.0; // under_relaxation factor
	
	// Testing
	//printf("v velocity from last iteration before updating coefficients \n");
	//print_v_velocity(M, N, velocity_star)
	
	// V velocity
	for(int i=0;i<M;i++) {
		for(int j=1;j<=N-1;j++) {
			// neighbouring coefficients parameters
					// Convective parameters
					Sp=0;
					Su=0;
					
					if(i==M-1) { // Right Boundary
						Fe=( BCs[0][2])*dy;
					}
					else {
						Fe=( velocity_star_old[i+1][j][0] + velocity_star_old[i+1][j-1][0] )*dy / 2.0;
					}
					if (i==0) { // Left Boundary
						Fw=( BCs[0][3] )*dy;
					}
					else {
						Fw=( velocity_star_old[i][j][0] + velocity_star_old[i][j-1][0] )*dy / 2.0;
					}
					Fn=( velocity_star_old[i][j+1][1] + velocity_star_old[i][j][1] )*dx / 2.0;
					Fs=( velocity_star_old[i][j][1] + velocity_star_old[i][j-1][1] )*dx / 2.0;
					// Diffusive parameters
					if (i==M-1) { // Right Boundary
						Sp = Sp - 2.0*dy / ( Re * dx);
						Su = Su + (2.0*BCs[1][2]*dy) / ( Re * dx);
					}
					else {
						De= 4*dy / ( 4 * Re * dx);
					}
					if (i==0) { // Left Boundary
						Sp = Sp - 2.0*dy / ( Re * dx);
						Su = Su + (2.0*BCs[1][3]*dy) / ( Re * dx);
					}
					else {
						Dw= 4*dy / ( 4 * Re * dx);
					}
					Dn= dx / ( Re * dy);
					Ds= dx / ( Re * dy);
					// neighbouring coefficients
					if (i==M-1) { // Right Boundary
						ae = 0;
					}
					else {
						ae = max(-Fe,(De-Fe/2.0),0);
					}
					if (i==0) { // Left Boundary
						aw = 0;
					}
					else {
						aw = max(Fw,(Dw+Fw/2.0),0);
					}
					an = max(-Fn,(Dn-Fn/2.0),0);
					as = max(Fs,(Ds+Fs/2.0),0);
					// Point term
					velocity_coefficients[i][j][1] = ae + aw + an + as + ( Fe - Fw + Fn - Fs  ) - Sp;
					ap = ae + aw + an + as + ( Fe - Fw + Fn - Fs  ) - Sp;
					// Pressure term coefficients or source term
					source = dx*(pressure_star[i][j-1] - pressure_star[i][j]) + Su;
					
					//vel = (ae*velocity_star[i+1][j][1] + aw*velocity_star[i-1][j][1]  +                              
					//      an*velocity_star[i][j+1][1] + as*velocity_star[i][j-1][1]  +    source ) / ap;
					
					// Solved using Gauss Seidel method
					vel = source;
					if(ae!=0) {vel = vel + ae*velocity_star[i+1][j][1];} 
					if(aw!=0) {vel = vel + aw*velocity_star[i-1][j][1];}
					if(an!=0) {vel = vel + an*velocity_star[i][j+1][1];}
					if(as!=0) {vel = vel + as*velocity_star[i][j-1][1];}   
					vel = vel / ap;
					velocity_star[i][j][1] =  (1-w)*velocity_star_old[i][j][1] + w*vel; 
					
					// Testing
					/*printf("v_velocity coefficients");
					printf("(%d,%d) --> ",i,j);
					printf("Fe=%lf , Fw=%lf , Fn=%lf , Fs=%lf   \n",Fe, Fw,Fn,Fs);
					printf("ap=%lf , ae=%lf , aw=%lf , an=%lf , as=%lf ,  source=%lf \n",velocity_coefficients[i][j][1][0], ae,aw,an,as,velocity_coefficients[i][j][1][5]);*/
		}
	}
	return 1;
}


double pressure_dash_update( int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity_star[M+1][N+1][2], double velocity_coefficients[M][N][2], double pressure_dash[M][N]  ) {
	double dx= grid_size[0];
	double dy= grid_size[1];
	
	// Initializing pressure equation coefficients
	double Fe=0, Fw=0, Fn=0, Fs=0; // Fe= east face, Fw= west face, Fn= north face, Fs= south face
	double de=0, dw=0, dn=0, ds=0; // De= east face, Dw= west face, Dn= north face, Ds= south face
	double ap=0, ae=0, aw=0, an=0, as=0, source=0 ; // ap=point value, ae= east face, aw= west face, an= north face, as= south face, source = source value
	
	double pres=0;
	double w=1.0; // Under_relaxation factor
	double error = 0;
	
	for(int i=0;i<M;i++) {
		for(int j=0;j<N;j++) {
				
				// pressure
					// neighbouring coefficients parameters
					// Convective parameters
					Fe=( velocity_star[i+1][j][0] )*dy;
					Fw=( velocity_star[i][j][0] )*dy;
					Fn=( velocity_star[i][j+1][1] )*dx;
					Fs=( velocity_star[i][j][1] )*dx;
					// d_IJ parameters
					de= dy / velocity_coefficients[i+1][j][0]; // From u velocity coeficient
					dw= dy / velocity_coefficients[i][j][0];
					dn= dx / velocity_coefficients[i][j+1][1]; // From v velocity coeficient
					ds= dx / velocity_coefficients[i][j][1];
					// neighbouring coefficients
					// east neighbour = ae
					if (i==M-1) { // East or right Boundary
						ae = 0;
					}
					else {
						ae = (dy) * de;
					}
					// west neighbour = aw
					if(i==0) { // West or Left Boundary
						aw = 0;
					}
					else {
						aw = (dy) * dw;
					}
					//  north neighbour = an
					if (j==N-1) { // North or Top Boundary
						an = 0;
					}
					else {
						an = (dx) * dn;
					}
					//  south neighbour = as
					if (j==0) { // South or Bottom Boundary
						as = 0;
					}
					else {
						as = (dx) * ds;
					}
					
					// Point term
					ap = ae + aw + an + as  ;
					// Pressure term coefficients or source term
					source = - ( Fe - Fw + Fn - Fs);
					error = error + pow(source,2);
				
					// Solved using Gauss Seidel method
					pres = source;
					if (ae!=0) {pres = pres + ae*pressure_dash[i+1][j];} 
					if (aw!=0) {pres = pres + aw*pressure_dash[i-1][j];}                                      	 
					if (an!=0) {pres = pres + an*pressure_dash[i][j+1];} 
					if (as!=0) {pres = pres + as*pressure_dash[i][j-1];} 
					pres = pres / ap;
					//printf("%lf ->",pres);
					pressure_dash[i][j] = (1-w)*pressure_dash[i][j] + w*pres;
					
					//printf("((i,j) => %lf) --\t",pressure_dash[i][j]);
				
					// Testing
					//printf("----Pressure Correction Equation Coefficients----\n");
					//printf("(%d,%d) --> ",i,j);
					//printf("ap = %lf \t",ap);
					//printf("%lf -->",pressure_dash[i][j]);
					//printf("Fe=%lf , Fw=%lf , Fn=%lf , Fs=%lf \n",Fe, Fw,Fn,Fs);
					//printf("ap=%lf , ae=%lf , aw=%lf , an=%lf , as=%lf , source=%lf \n",ap, ae,aw,an,as,source);

		}
	}
	error = pow( ( error / (M*N) ) , 0.5 ); 
	return error;
}



double max(double a, double b, double c) {
	double result = a;
	if (b>result) {
		result=b;
	}
	if (c>result) {
		result=c;
	}
	return result;
}

double Results_Collocated_Grid(int M, int N, double Re, double grid_size[2], double BCs[2][4], double velocity[M+1][N+1][2]) {
	// Shifting to collocated Grid
	double u_velocity[M+1][N+1];
	double v_velocity[M+1][N+1];
	double StreamFunction[M+1][N+1];
	double vorticity[M+1][N+1];
	
	double dx= grid_size[0];
	double dy= grid_size[1];
	
	printf("-------- Calculating All Results on Collocated Grid --------\n");
	
	for(int i=0;i<=M;i++) {
		for (int j=0;j<=N;j++) { // For Corner points Boundary
			if ((i==0 && j==0) || (i==0 && j==N) || (i==M && j==0) || (i==M && j==N)) {
				continue; // Corner points
			}
			else if(j==0) { // For Bottom Boundary
				u_velocity[i][j] = BCs[0][1];
				v_velocity[i][j] = BCs[1][1];
			}
			else if(j==N) { // For Top Boundary
				u_velocity[i][j] = BCs[0][0];
				v_velocity[i][j] = BCs[1][0];
			}
			else if(i==0) { // For Left Boundary
				u_velocity[i][j] = BCs[0][3];
				v_velocity[i][j] = BCs[1][3];
			}
			
			else if(i==M) { // For Right Boundary
				u_velocity[i][j] = BCs[0][2];
				v_velocity[i][j] = BCs[1][2];
			}
			else {
				u_velocity[i][j] = ( velocity[i][j][0] +  velocity[i][j-1][0] ) / 2.0;
				v_velocity[i][j] = ( velocity[i][j][1] +  velocity[i-1][j][1] ) / 2.0;
			}
		}
	}
	// Corner_points
	u_velocity[0][0] = ( u_velocity[0][1] + u_velocity[1][0] ) / 2.0;  	v_velocity[0][0] = ( v_velocity[0][1] + v_velocity[1][0] ) / 2.0;
	u_velocity[0][N] = ( u_velocity[0][N-1] + u_velocity[1][N] ) / 2.0;	v_velocity[0][N] = ( v_velocity[0][N-1] + v_velocity[1][N] ) / 2.0;
	u_velocity[M][0] = ( u_velocity[M-1][0] + u_velocity[M][1] ) / 2.0;	v_velocity[M][0] = ( v_velocity[M-1][0] + v_velocity[M][1] ) / 2.0;
	u_velocity[M][N] = ( u_velocity[M-1][N] + u_velocity[M][N-1] ) / 2.0;	v_velocity[M][N] = ( v_velocity[M-1][N] + v_velocity[M][N-1] ) / 2.0;
	
	printf("-------- u and v velocities calculated on Collocated Grid --------\n");
	
	// Vorticity Calculation
	for(int i=0;i<=M;i++) {
		for(int j=0;j<=N;j++){
			vorticity[i][j] = 0;
			if(i==0) { // For Left Boundary
				// vorticity[i][j] = vorticity[i][j] + ( - v_velocity[i+2][j] + 4*v_velocity[i+1][j] - 3*v_velocity[i][j] ) / (2.0*dx);
				vorticity[i][j] = vorticity[i][j] + (v_velocity[i+1][j] - v_velocity[i][j] ) / (dx);
			}
			if(i==M) { // For Right Boundary
				//vorticity[i][j] = vorticity[i][j] + ( v_velocity[i-2][j] - 4*v_velocity[i-1][j] + 3*v_velocity[i][j] ) / (2.0*dx);
				vorticity[i][j] = vorticity[i][j] + (v_velocity[i][j] - v_velocity[i-1][j] ) / (dx);
			}
			if(j==0) { // For Bottom Boundary
				//vorticity[i][j] = vorticity[i][j] -1*( - u_velocity[i][j+2] + 4*u_velocity[i][j+1] - 3*u_velocity[i][j] ) / (2.0*dy);
				vorticity[i][j] = vorticity[i][j] - ( u_velocity[i][j+1] - u_velocity[i][j] ) / (dy);
			}
			if(j==N) { // For Top Boundary
				//vorticity[i][j] = vorticity[i][j] -1*( u_velocity[i][j-2] - 4*u_velocity[i][j-1] + 3*u_velocity[i][j] ) / (2.0*dy);
				vorticity[i][j] = vorticity[i][j] - ( u_velocity[i][j] - u_velocity[i][j-1] ) / (dy);
			}
			if (i!=0 && j!=0 && i!=M && j!=N) {
				vorticity[i][j] = ( (v_velocity[i+1][j] - v_velocity[i-1][j]) / (2.0*dx) ) - ( (u_velocity[i][j+1] - u_velocity[i][j-1]) / (2.0*dy) );
			}
		}
	}
	printf("-------- Vorticity calculated on Collocated Grid --------\n");
	
	// StreamFunction Calculation
	double beta = dx/dy;
	double epsilon = pow(10,-6);
	double psi_error=1, Sf=0;
	double psi_min=0;
	int psi_min_i=0, psi_min_j=0;
	int iter_count=0;
	// Defining Boundaries
	for(int j=0;j<=N;j++) {
		StreamFunction[0][j] = 0;
		StreamFunction[M][j] = 0;
	}
	for(int i=0;i<=M;i++){
		StreamFunction[i][0] = 0;
		StreamFunction[i][N] = 0;
	}
	printf("-------- Calculating Stream Function on Collocated Grid using Gauss Seidel Method --------\n");
	while (psi_error>epsilon) {
		psi_error=0; // L2 norm
		for(int i=1;i<M;i++) {
			for(int j=1;j<N;j++){ // Gauss Seidel Method
				Sf = ( 0.5/(1.0+pow(beta,2)) ) * ( vorticity[i][j]*pow(dx,2) + pow(beta,2)*(StreamFunction[i][j+1] + StreamFunction[i][j-1]) + 			                                                                                     StreamFunction[i+1][j] + StreamFunction[i-1][j] );
				psi_error = psi_error + pow(StreamFunction[i][j]-Sf,2);
				StreamFunction[i][j] = Sf;
			}
		}
		psi_error=pow(psi_error/((M-1)*(N-1)),0.5);
		iter_count++;
	}
	printf("The Iteration for Calculating Stream Function are as follows: %d \n",iter_count);
	printf("-------- Stream Dunction calculated on Collocated Grid --------\n");
	
	// Finding minimum streamfunction value and location
	for(int i=1;i<M;i++) {
			for(int j=1;j<N;j++){ 
				if (psi_min>StreamFunction[i][j]) {
					psi_min = StreamFunction[i][j];
					psi_min_i = i;
					psi_min_j = j;
				}
			}
	}
	printf("The minimum stream function value is %lf located at (%lf, %lf) \n", psi_min,psi_min_i*dx,psi_min_j*dy);
	printf("The Vorticity value is %lf at (%lf, %lf) \n", vorticity[psi_min_i][psi_min_j],psi_min_i*dx,psi_min_j*dy);
	
	// printing Results
	double x,y;
	char str[10];
	sprintf(str, "%d", (int)Re);
	char file_name[100] = "Results_Re_";
	strcat(file_name,str);
	strcat(file_name,".dat");
	FILE *fp1, *fp6;
	fp1 = fopen(file_name,"w");
	fprintf(fp1,"ZONE I=%d, J=%d \n",M,N);
	for(int j=0;j<=N;j++) {
		y=(j)*(dy);
		for(int i=0;i<=M;i++) {
			x=(i)*dx;
			fprintf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u_velocity[i][j],v_velocity[i][j],StreamFunction[i][j],vorticity[i][j]);
		}
	}
		
	
	// Generating Results
	FILE *fp2,*fp3,*fp4,*fp5;
	char fname_1[100] = "u_velocity_centerline_x_Re_";
	char fname_2[100] = "u_velocity_centerline_y__Re_";
	char fname_3[100] = "v_velocity_centerline_x_Re_";
	char fname_4[100] = "v_velocity_centerline_y_Re_";
	strcat(fname_1,str);	strcat(fname_1,".dat");
	strcat(fname_2,str);	strcat(fname_2,".dat");
	strcat(fname_3,str);	strcat(fname_3,".dat");
	strcat(fname_4,str);	strcat(fname_4,".dat");
		
	fp2 = fopen(fname_1,"w");
	fp3 = fopen(fname_2,"w");
	fp4 = fopen(fname_3,"w");
	fp5 = fopen(fname_4,"w");
	// u_velocity
	for(int j=0;j<=N;j++) { //x=L/2
		fprintf(fp2,"%lf %lf\n",u_velocity[64][j],j*dy);
	}
	for(int i=0;i<=M;i++) { //y=L/2
		fprintf(fp3,"%lf %lf\n",i*dx, u_velocity[i][64]);
	}
	// v_velocity
	for(int j=0;j<=N;j++) { //x=L/2
		fprintf(fp4,"%lf %lf\n",v_velocity[64][j],j*dy);
	}
	for(int i=0;i<=M;i++) { //y=L/2
		fprintf(fp5,"%lf %lf\n",i*dx, v_velocity[i][64]);
	}
			
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	
	return 1;
}


void print_u_velocity(int M, int N, double vel[M+1][N+1][2]) {
	printf("\n");
	printf("----------xxxxxxxxxxxxx-------U Velocity-------xxxxxxxxxxxxxxx--------------\n\n");
	for (int j=N-1;j>=0;j--) {
		printf("%d --> \n",j);
		for(int i=0;i<=M;i++) {
			printf("%lf \t",vel[i][j][0]);
		}
		printf("\n");
	}
	printf("----------xxxxxxxxxxxxx--------------xxxxxxxxxxxxxxx--------------\n\n");
}

void print_v_velocity(int M, int N, double vel[M+1][N+1][2]) {
	printf("\n");
	printf("----------xxxxxxxxxxxxx-------V Velocity-------xxxxxxxxxxxxxxx--------------\n\n");
	for (int j=N;j>=0;j--) {
		printf("%d --> \n",j);
		for(int i=0;i<M;i++) {
			printf("%lf \t",vel[i][j][1]);
		}
		printf("\n");
	}
	printf("----------xxxxxxxxxxxxx--------------xxxxxxxxxxxxxxx--------------\n\n");
}

void print_pressure(int M, int N, double pres[M][N]) {
	printf("\n");
	printf("----------xxxxxxxxxxxxx-------Pressure-------xxxxxxxxxxxxxxx--------------\n\n");
	for (int j=N-1;j>=0;j--) {
		printf("%d --> \n",j);
		for(int i=0;i<M;i++) {
			printf("%lf \t",pres[i][j]);
		}
		printf("\n");
	}
	printf("----------xxxxxxxxxxxxx--------------xxxxxxxxxxxxxxx--------------\n\n");
}
