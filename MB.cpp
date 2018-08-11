#include <iostream>
#include <fstream>
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

std::ofstream myfile;


/************************************* FUNCTIONS ***************************************************/

//clear old data from output file
void wipeOutput(){
	myfile.open("output.txt", std::ofstream::out | std::ofstream::trunc);
	myfile.close();

};

void writeConstantsToFile(int N, int T, int p, double dt, double dx){
	
	myfile.open ("output.txt", std::fstream::app);
	myfile << N << " " << T << " " << p << " " << dt << " " << dx;
	myfile.close();
}
/*
enforces periodic boundary conditions for the discreet laplacian
*/
int BC(int i){
	
	int I = (i + N) % N;

	if (I > N - 1)
		cout << "ERROR!!";

	if (I < 0)
		cout << "ERROR!!";

	return I;

};
//gets values in the phi array
double phiVal(int i, int j){
	ptrPHI = begining_Of_PHI;
	double element = * (ptrPHI + (BC(i) * N) + BC(j));
	return element;

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


//evaluates the laplacian for phi(i,j). Eqn B.5 in Nik's book may be used as a more isotropic laplacian
//for PERIODIC BOUNDARY CONDITIONS - right now doesn't deal with 0 nodes well...
double Laplacian_phi(int i, int j){

	double laplacian_ij; 


	laplacian_ij = phiVal(i + 1 , j) + phiVal( i - 1, j)
			+ phiVal( i, i + 1 ) + phiVal(i, i - 1) 
			- 4 * ( phiVal(i,j) );
	
	return laplacian_ij;
};

double MU(int i, int j){
	double mu = - ( W / pow(dx, 2) ) * Laplacian_phi(i,j) + df(i,j);
}

double Laplacian_mu(int i, int j){
	double newMU;
	newMU = MU(i + 1, j) + MU(i - 1, j) + MU(i, j - 1) + MU(i, j + 1) - 4 * MU(i,j);
	return newMU;
};


//single time march for a type A model
//only stable for sufficiently small time steps (delta_t)
double timeMarchA(int i, int j){

	double newPHI = phiVal(i,j) + (dt_bar / pow(dx_bar, 2)) * Laplacian_phi(i , j)
						- dt_bar * df(i,j);

	return newPHI;

};

double timeMarchB(int i, int j){
	double mu, newPHI;

	newPHI = phiVal(i,j) + ( dt * M / pow(dx, 2) ) * Laplacian_mu(i, j);


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
	dt = 0.05; 
	dx = 0.8;

	//assign vaue for N
	N = 100;
	T = 100; //make sure T is a multiple of 10 to make timesteps are good.
	int skipPrint = 10;

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
		if( timestep % skipPrint == 0 ){
			printPHI();
		}

	} 

	//stop stopwatch
	auto finish = std::chrono::high_resolution_clock::now();

	//should print parameters and other info to top of the output file.
	writeConstantsToFile(N, T, skipPrint, dt, dx);

	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
        std::cout << microseconds.count() << "µs\n";
        if (microseconds > std::chrono::seconds(1))

	return 0;
}