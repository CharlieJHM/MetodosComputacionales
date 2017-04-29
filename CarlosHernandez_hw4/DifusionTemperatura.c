#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//Declaracion de funciones -----------------------------------------------------------------
void ini_condition(double *T, int n_points, double dx, int caso, int boundary);
void copy(double *T_to_copy, double *T_copy, int n_points);
void move_dt_1(double *T_now, double *T_bef, int n_points, double alpha, int bound, int caso);
double mean_T(double *T_to_prom, int n_points);
void print_matrix(double *T, int n_points);
void append_matrix_to_file(FILE in, double *Matrix, int n_points);

//##########################################################################################
//Main -------------------------------------------------------------------------------------

int main(){
	double x_length = 1.0;
	double y_lenght = 1.0;
	double dx = 0.01;
	int n_x = x_length/dx;

	double t_max = 2500.0;
	double dt = 0.1;
	int n_t = t_max/dt;

	double nu = 1E-4;
	double alpha = nu*dt/(dx*dx);

	FILE *in;
	FILE *in2;
	char filename[100] = "Tmean_time.txt";
	char filename2[100] = "T_matrices_anim.txt";

	double *T_past = malloc((n_x*n_x)*sizeof(double));
	double *T_present = malloc((n_x*n_x)*sizeof(double));
	double *T_mean = malloc(n_t*sizeof(double));

	//Caso 1, Boundary abierto -----------------------------------------------
	ini_condition(T_past, n_x, dx, 1, 1);
	int i;
	double t;
	for(i=0;i<n_t;i++){
		t = i*dt;
		if(t == 0.0 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 1, 1);
		copy(T_present, T_past, n_x);
	}

	in = fopen(filename, "w");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Caso 1, Boundary periodico --------------------------------------------
	ini_condition(T_past, n_x, dx, 1, 2);
	in2 = fopen(filename2, "w");
	for(i=0;i<n_t;i++){
		t = i*dt;
		if((i%20==0 && i<=460) || (i%500==0 && 500<=i && i<=5500)){
			int k;
			int j = 0;
			for(k=0;k<(n_x*n_x);k++){
				if(j == 99){
					fprintf(in2, "%f\n", T_past[k]);
					j = 0;
				}else{
					fprintf(in2, "%f\t", T_past[k]);
					j = j + 1;
				}
			}
		}

		if(t == 0.0 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 2, 1);
		copy(T_present, T_past, n_x);
	}
	fclose(in2);

	in = fopen(filename, "a");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Caso 1, Boundary fijo----------------------------------------------------
	ini_condition(T_past, n_x, dx, 1, 3);
	for(i=0;i<n_t;i++){
		t = i*dt;
		if(t == 0.0 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 3, 1);
		copy(T_present, T_past, n_x);
	}

	in = fopen(filename, "a");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Caso 2, Boundary abierto -----------------------------------------------
	ini_condition(T_past, n_x, dx, 2, 1);
	for(i=0;i<n_t;i++){
		t = i*dt;
		if(t == 0.0 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 1, 2);
		copy(T_present, T_past, n_x);
	}

	in = fopen(filename, "a");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Caso 2, Boundary periodico --------------------------------------------
	ini_condition(T_past, n_x, dx, 2, 2);
	for(i=0;i<n_t;i++){
		t = i*dt;
		if(t == 0.0 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 2, 2);
		copy(T_present, T_past, n_x);
	}

	in = fopen(filename, "a");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Caso 2, Boundary fijo----------------------------------------------------
	ini_condition(T_past, n_x, dx, 2, 3);
	for(i=0;i<n_t;i++){
		t = i*dt;
		if(t == 0.1 || t == 100.0 || t == 2499.9){
			print_matrix(T_past, n_x);
			printf("\n");
		}
		T_mean[i] = mean_T(T_past, n_x);
		move_dt_1(T_present, T_past, n_x, alpha, 3, 2);
		copy(T_present, T_past, n_x);
	}

	in = fopen(filename, "a");
	for(i=0;i<n_t;i++){
		if(i == n_t-1){
			fprintf(in, "%f\n", T_mean[i]);
		}else{
			fprintf(in, "%f\t", T_mean[i]);
		}
	}
	fclose(in);

	//Liberacion de memoria----------------------------------------------------
	free(T_past);
	free(T_present);
	free(T_mean);
}

//##########################################################################################
//Funcion de inicializacion ------------------------------------------------------
//caso = 1, 2 corresponden al caso 1, caso 2 de la tarea respectivamente
//boundary = 1, 2, 3 correponden a condiciones de frontera abiertas, periodicas y fijas respectivamente

void ini_condition(double *T, int n_points, double dx, int caso, int boundary){
	int i;
	int j = 0;
	int k = 0;
	if(caso == 1 || caso == 2){
		for(i=0;i<(n_points*n_points);i++){
			double x = j*dx;
			double y = k*dx;
			if(0.3 <= x && x <= 0.5 && 0.45 <= y && y <= 0.55){
				T[i] = 100.0;
				k = k + 1;
			}else{
				T[i] = 50.0;
				k = k + 1;
			}

			if(k==100){
				j = j + 1;
				k = 0;
			}
		}
	}
}

//##########################################################################################
//Funcion que copia los valores de un array a otro -------------------------------

void copy(double *T_to_copy, double *T_copy, int n_points){
	int i;
	for(i=0;i<(n_points*n_points);i++){
		T_copy[i] = T_to_copy[i];
	}
}

//##########################################################################################
//Funcion que evoluciona el sistema en un dt para --------------------------------------------------
//bound = 1, 2, 3 correponden a condiciones de frontera abiertas, periodicas y fijas respectivamente

void move_dt_1(double *T_now, double *T_bef, int n_points, double alpha, int bound, int caso){
	int i;
	int j = 0;
	int k = 0;
	if(bound == 1){
		for(i=0;i<(n_points*n_points);i++){
			if(j==0 && k==0){
				T_now[i] = alpha*(T_bef[1]+T_bef[100]) + (1.0 - 4.0*alpha)*T_bef[0];
				k = k + 1;
			}else if(j==0 && k==99){
				T_now[i] = alpha*(T_bef[98]+T_bef[199]) + (1.0 - 4.0*alpha)*T_bef[99];
				k = k + 1;
			}else if(j==99 && k==0){
				T_now[i] = alpha*(T_bef[9800]+T_bef[9901]) + (1.0 - 4.0*alpha)*T_bef[9900];
				k = k + 1;
			}else if(j==99 && k==99){
				T_now[i] = alpha*(T_bef[9899]+T_bef[9998]) + (1.0 - 4.0*alpha)*T_bef[9999];
				k = k + 1;
			}else if(k==0){
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i+1]+T_bef[i+100]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(j==0){
				T_now[i] = alpha*(T_bef[i-1]+T_bef[i+100]+T_bef[i+1]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(k==99){
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i-1]+T_bef[i+100]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(j==99){
				T_now[i] = alpha*(T_bef[i-1]+T_bef[i-100]+T_bef[i+1]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(caso==2 && 30<=j && j<=50 && 45<=k && k<=55){
				T_now[i] = 100.0;
				k = k + 1;
			}else{
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i+1]+T_bef[i+100]+T_bef[i-1]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}
			if(k==100){
				j = j + 1;
				k = 0;
			}
		}
	}

	if(bound == 2){
		for(i=0;i<(n_points*n_points);i++){
			if(j==0 && k==0){
				T_now[i] = alpha*(T_bef[1]+T_bef[100]+T_bef[99]+T_bef[9900]) + (1.0 - 4.0*alpha)*T_bef[0];
				k = k + 1;
			}else if(j==0 && k==99){
				T_now[i] = alpha*(T_bef[98]+T_bef[199]+T_bef[0]+T_bef[9999]) + (1.0 - 4.0*alpha)*T_bef[99];
				k = k + 1;
			}else if(j==99 && k==0){
				T_now[i] = alpha*(T_bef[9800]+T_bef[9901]+T_bef[0]+T_bef[9999]) + (1.0 - 4.0*alpha)*T_bef[9900];
				k = k + 1;
			}else if(j==99 && k==99){
				T_now[i] = alpha*(T_bef[9899]+T_bef[9998]+T_bef[9900]+T_bef[99]) + (1.0 - 4.0*alpha)*T_bef[9999];
				k = k + 1;
			}else if(k==0){
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i+1]+T_bef[i+100]+T_bef[i+99]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(j==0){
				T_now[i] = alpha*(T_bef[i-1]+T_bef[i+100]+T_bef[i+1]+T_bef[i+9900]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(k==99){
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i-1]+T_bef[i+100]+T_bef[i-99]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(j==99){
				T_now[i] = alpha*(T_bef[i-1]+T_bef[i-100]+T_bef[i+1]+T_bef[i-9900]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}else if(caso==2 && 30<=j && j<=50 && 45<=k && k<=55){
				T_now[i] = 100.0;
				k = k + 1;
			}else{
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i+1]+T_bef[i+100]+T_bef[i-1]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}
			if(k==100){
				j = j + 1;
				k = 0;
			}
		}
	}

	if(bound == 3){
		for(i=0;i<(n_points*n_points);i++){
			if(j==0 || k==0 || j==99){
				T_now[i] = 50.0;
				k = k + 1;
			}else if(caso==2 && 30<=j && j<=50 && 45<=k && k<=55){
				T_now[i] = 100.0;
				k = k + 1;
			}else{
				T_now[i] = alpha*(T_bef[i-100]+T_bef[i+1]+T_bef[i+100]+T_bef[i-1]) + (1.0 - 4.0*alpha)*T_bef[i];
				k = k + 1;
			}
			if(k==100){
				j = j + 1;
				k = 0;
			}
		}
	}
}

//##########################################################################################
//Funcion que halla el valor promedio del array que se le entrega ----------------

double mean_T(double *T_to_prom, int n_points){
	int i;
	double mean;
	double sum = 0.0;
	for(i=0;i<(n_points*n_points);i++){
		sum = sum + T_to_prom[i];
	}
	mean = sum/(n_points*n_points);
	return mean;
}

//##########################################################################################
//Funcion que imprime en forma de matriz de n_points*n_points el vector que se le entregue
void print_matrix(double *T, int n_points){
	int i;
	int j = 0;
	for(i=0;i<(n_points*n_points);i++){
		if(j == 99){
			printf("%f\n", T[i]);
			j = 0;
		}else{
			printf("%f\t", T[i]);
			j = j + 1;
		}
	}
}

//##########################################################################################
//Funcion que fija la temperatura del rectangulo inicial a 100C
void caso2(double *T, int n_points, double dx){
	int i;
	int j = 0;
	int k = 0;
	for(i=0;i<(n_points*n_points);i++){
		double x = j*dx;
		double y = k*dx;
		if(0.3 <= x && x <= 0.5 && 0.45 <= y && y <= 0.55){
			T[i] = 100.0;
			k = k + 1;
		}

		if(k==100){
			j = j + 1;
			k = 0;
		}
	}
}
