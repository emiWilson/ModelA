#include <iostream>
#include <fstream>
#include <functional>
#include <math.h> /*pow*/
#include <random>
#include <cstdlib>
#include <chrono> //time of execution of program, want to scale with system size t^(1/2)
 //PROBLEM MAY BE IN INDEXING.

using namespace std;

//is there a way to do this without hardcoding the phi array size??????

/************************************ GLOBAL VARIABLES ******************************************/

// N is the number of grid points in the numerical array. N = L / (delta_x), where L is the domain length in each direction
int N, T; //number of physical and temporal iterations


//ptrPHI will be used to navigate the NxN phi array: phi[i][j] = *(ptr + (i * N) + j))
double * ptrPHI ;

//this pointer will always point to the start of the PHI array and will remain unaltered. Can't make "const" variable type b/c don't know the phi array before the program starts. 
double * begining_Of_PHI;

double dx, dy, dt;

double dx_bar, dt_bar;

double W, M;

int above_ij, below_ij, left_ij, right_ij; //stand for above, below, left, right

std::ofstream myfile;


/************************************* FUNCTIONS ***************************************************/

//clear old data from output file
void wipeOutput(){
	myfile.open("output.txt", std::ofstream::out | std::ofstream::trunc);
	myfile.close();

};

void writeConstantsToFile(int T, int N, int T_print, double dt, double dx){
	
	myfile.open ("output.txt", std::fstream::app);
	myfile << T << " " << N << " " << T_print << " " << dt << dx;
	myfile.close();
}

//gets values in the phi array
double phiVal(int i, int j){
	ptrPHI = begining_Of_PHI;
	double element = * (ptrPHI + (i * N) + j);
	return element;

};

void periodic_BC(int i, int j){
	above_ij = i + 1;
	below_ij = i - 1;
	left_ij = j + 1;
	right_ij = j - 1;

	//cout << "laplacian loop" << phiVal(0,0);

	//periodic boundary conditions
	if (i == 0){
		below_ij = 0;
	}
	if (i == N - 1){
		above_ij = 0;
	}
	if (j == 0){
		right_ij = 0;
	}
	if (j == N - 1){
		left_ij = 0;
	}
};


//evaluates the laplacian for phi(i,j). Eqn B.5 in Nik's book may be used as a more isotropic laplacian
//for PERIODIC BOUNDARY CONDITIONS - right now doesn't deal with 0 nodes well...
double Laplacian(int i, int j){

	double laplacian_ij; 

	periodic_BC(i,j);


	laplacian_ij = phiVal( above_ij, j) + phiVal( below_ij, j)
			+ phiVal( i, left_ij ) + phiVal(i, right_ij) 
			- 4 * ( phiVal(i,j) );
	//cout << " " << phiVal( above, j) << " " << phiVal(below, j)
	//		<< " " << phiVal( i, left ) << " " << phiVal(i, right) << " ij is " << i << " " << j << endl;
	
	return laplacian_ij;
};


/**
function for the derivative of the the free energy with respect to the order parameter
df/dPhi

f(phi) will give a double well potential where the two minima correspond to the solid and liquid phases

should plot f(phi) with the values for the constants to see what it looks like.
**/

double df(int i, int j){

	double a2 = -1;
	double a4 = 1;


	double diff_f = a2 * phiVal(i, j) + a4 * pow( phiVal(i,j) , 3) ; // + higher order terms - EQ 2.38 in book.

	return diff_f;
};

//single time march for a type A model
//only stable for sufficiently small time steps (delta_t)
double timeMarchA(int i, int j){

	double newPHI = phiVal(i,j) + (dt_bar / pow(dx_bar, 2)) * Laplacian(i , j)
						- dt_bar * df(i,j);

	//cout << newPHI << " from " << phiVal(i,j) << endl;
	
	return newPHI;

};


double Laplacian_mu(int i, int j){
	double mu_ij, laplacian_phi, laplacian_f; 

	periodic_BC(i,j);

	laplacian_phi = Laplacian(above_ij, j) + Laplacian(below_ij, j)
			+ Laplacian(i, left_ij ) + Laplacian(i, right_ij) 
			- 4 * (Laplacian(i,j));
	
	laplacian_f = df(above_ij, j) + df(below_ij, j)
			+ df(i, left_ij) + df(i, right_ij)
			-4 * (df(i,j));

	mu_ij = - laplacian_phi / ( pow(dx_bar, 2) ) + laplacian_f;

	return mu_ij;
	
};

double timeMarchB(int i, int j){

	double newPHI = phiVal(i,j) + (dt_bar / pow(dx_bar,2)) * Laplacian_mu(i,j);
	

	return newPHI;
};





//be careful not to print this array every time step. it will quickly fill up disk space
void printPHI(){

  	myfile.open ("output.txt", std::fstream::app);

	for (int i = 0; i < N; i ++){
		for (int j = 0; j < N; j ++){
			myfile << phiVal(i,j) << " ";
		}
		myfile << endl;
	}

	myfile << endl;

	myfile.close();

};





/********************************* MAIN METHOD *****************************************************/
int main(){

	wipeOutput(); //ensure all old data is deleted from file.

	//magnitude of timestep. should obey stability restriction delta_t <  (delta_x ^ 2) / 4
	//restriction means that it is not possible to advance a solution explicity faster than the inherent diffusion
	//time of the problem.
	dt = 0.1; 
	dx = 0.8;

	//assign vaue for N
	N = 100;
	T = 10; //make sure T is a multiple of 10 to make timesteps are good.
	int T_print = 1;

	//other variables
	W = pow(0.25, 1/2);
	M = 1; //lookup realistic value for M.
	dt_bar = dt * M;
	dx_bar = dx / W;


	//start stopwatch
	auto start = std::chrono::high_resolution_clock::now();

	/********************************************
			create phi array
	*********************************************/
	double phi[N][N];

	begining_Of_PHI = & phi[0][0]; //point ptrPHI to start of 2D phi array
	ptrPHI = begining_Of_PHI;



	/********************************************
			Fill phi array
	********************************************/

	std::default_random_engine de(time(0));
  	std::normal_distribution<double> distribution(0, 0.001); // 0 mean and 0.001 standard deviation
	
	for (int i = 0; i < N ; i ++){
		for (int j = 0; j < N; j ++){
			
			phi[i][j] = distribution(de);
	
		}
	}
	

	printPHI();


	
	for(int timestep = 0; timestep < T; timestep++){

		for (int i = 0; i < N; i ++){
			for (int j = 0; j < N ; j ++){
				phi[i][j] = timeMarchB(i,j);
			}
		}

		//print phi array for every 10th timestep
		if( timestep % T_print == 0 ){
			printPHI();
		}

	} 

	//stop stopwatch
	auto finish = std::chrono::high_resolution_clock::now();

	//should print parameters and other info to top of the output file.
	writeConstantsToFile(T, N, T/T_print, dt, dx);

	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
        std::cout << microseconds.count() << "µs\n";
        if (microseconds > std::chrono::seconds(1))

	return 0;
}