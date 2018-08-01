#include <iostream>
#include <fstream>
#include <math.h> /*pow*/
#include <cstdlib>

using namespace std;

//is there a way to do this without hardcoding the phi array size??????

/************************************ GLOBAL VARIABLES ******************************************/

// N is the number of grid points in the numerical array. N = L / (delta_x), where L is the domain length in each direction
int N;


//ptrPHI will be used to navigate the NxN phi array: phi[i][j] = *(ptr + (i * N) + j))
double * ptrPHI;

//this pointer will always point to the start of the PHI array and will remain unaltered. Can't make "const" variable type b/c don't know the phi array before the program starts. 
double * begining_Of_PHI;

double dx, dy, dt;

double T, Tm, H, L;


/************************************* FUNCTIONS ***************************************************/

//clear old data from output file
void wipeOutput(){
	std::ofstream ofs;
	ofs.open("output.txt", std::ofstream::out | std::ofstream::trunc);
	ofs.close();

};

//gets values in the phi array
double phiVal(int i, int j){
	ptrPHI = begining_Of_PHI;
	double element = * (ptrPHI + (i * N) + j);
	return element;

};



//evaluates the laplacian for phi(i,j). Eqn B.5 in Nik's book may be used as a more isotropic laplacian
//for PERIODIC BOUNDARY CONDITIONS - right now doesn't deal with 0 nodes well...
double Laplacian(int i, int j){

	double laplacian_ij; 

	int above = i + 1;
	int below = i - 1;
	int left = j + 1;
	int right = j - 1;

	//periodic boundary conditions
	if (i == 0){
		below = N;
	}
	if (i == N){
		above = 0;
	}
	if (j == 0){
		right = N;
	}
	if (j == N){
		left = 0;
	}

	laplacian_ij = phiVal( above, j) + phiVal( below, j)
			+ phiVal( i, left ) + phiVal(i, right) 
			- 4 * ( phiVal(i,j) );
	

	return laplacian_ij;
};

/**
function for the derivative of the the free energy with respect to the order parameter
df/dPhi

f(phi) will give a double well potential where the two minima correspond to the solid and liquid phases

should plot f(phi) with the values for the constants to see what it looks like.
**/
double df(int i, int j){

	double A = 2 * H + 6 * (T / Tm - 1) * L ;
	double B = 6 * H + (T / Tm - 1) * L;
	double C = 4 * H;

	double diff_f = A * phiVal(i, j ) + C * pow( phiVal(i,j) , 3);

	return diff_f;
};

//single time march.
//only stable for sufficiently small time steps (delta_t)
double timeMarch(int i, int j){

	double newPHI = phiVal(i,j) + (dt / pow(dx, 2)) * Laplacian(i , j)
						- dt * df(i,j);
	
	return newPHI;

};


//be careful not to print this array every time step. it will quickly fill up disk space
void printPHI(){
	ofstream myfile;
	myfile.open ("output.txt", std::fstream::app);

  	

	for (int i = 0; i < N; i ++){
		for (int j = 0; j < N; j ++){
			myfile << *(ptrPHI + (i + N) + j) << " ";
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
	dt = 0.05; // values given in Nik's 2007 paper
	dx = 1.1;
	dy = dx; //for simplicity, not required. 

	//assign value for N
	N = 20;

	//other variables
	T = 20; //the temperature, here considered a constant
	Tm = 50; //the melting temperature of the material
	H = 4; //the nucleation barrier
	L = 5.0; //latent heat of fusion


	//should print parameters and other info to top of the output file.

	/********************************************
			create phi array
	*********************************************/
	double phi[N][N];

	begining_Of_PHI = & phi[0][0]; //point ptrPHI to start of 2D phi array
	ptrPHI = begining_Of_PHI;


	/********************************************

	create time array, will be filled with integers from 0 ... N-1
	//same for x and y arrays

	********************************************/
	int t[N];
	int x[N], y[N];

	for (int i = 0; i < N; i ++){
		t[i] = i;
		x[i] = i;
		y[i] = i;
	}

	/********************************************
			Fill phi array
	********************************************/
	
	srand (time(NULL));
	double seed;

	for (int i = 0; i < N ; i ++){
		for (int j = 0; j < N; j ++){
			//phi[i][j] = (double) (i + 1) * (j + 1) ;
			seed = ((double)rand()) / ((double)RAND_MAX) * 1.0;
			//fill phi array with random 1s and 0s
			phi[i][j] = seed	;	
		}
	}

	printPHI();

	for(int timestep = 0; timestep < N * N; timestep++){

		for (int i = 0; i < N; i ++){
			for (int j = 0; j < N ; j ++){
				phi[i][j] = timeMarch(i,j);
			}
		}

		//print phi array for every 10th
		if( timestep % 10 == 0 ){
			printPHI();
		}

	}




	return 0;
}