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

double U, D_bar, a_s, epsilon_prime, a12, epsilon_4, lambda_bar; //Model C variables
double Der;

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

//explicit form of the x and y derivatives evaluated at right and left edges are
double DERX_iPlusHalf(int i, int j){

	Der = ( 1/dx_bar ) * (phiVal(i + 1, j) - phiVal(i,j));
	return Der;

};

double DERX_iMinusHalf(int i, int j){

	Der = (1 / dx_bar) * (phiVal(i,j) - phiVal(i - 1, j));
	return Der; 

};
double DERX_jPlusHalf(int i, int j){

	Der = ( 1/dx_bar ) * (phiVal(i, j + 1) - phiVal(i,j));
	return Der;

};

double DERX_jMinusHalf(int i, int j){

	Der = (1 / dx_bar) * (phiVal(i,j) - phiVal(i, j - 1));
	return Der; 

};


//"for the y derivative on the right edge of the finite volume, interpolation from the next and next
//nearest neighbout must be used"
double DERY_iPlusHalf(int i, int j){
	Der = (1 / ( 4 * dx_bar )) * 
		( phiVal(i + 1, j + 1) + phiVal(i, j + 1) + phiVal(i,j) + phiVal(i + 1, j) 
		- (phiVal(i + 1, j) + phiVal(i,j) + phiVal(i, j - 1) + phiVal(i + 1, j - 1)));
	return Der;

};

double DERY_iMinusHalf(int i, int j){
	Der = (1 / ( 4 * dx_bar )) * 
		( phiVal(i, j + 1) + phiVal(i - 1, j + 1) + phiVal(i - 1, j) + phiVal(i, j) 
		- (phiVal(i, j) + phiVal(i - 1,j) + phiVal(i - 1, j - 1) + phiVal(i, j - 1)));
	return Der; 
};

double DERY_jPlusHalf(int i, int j){
	Der = (1 / ( 4 * dx_bar )) * 
		( phiVal(i + 1, j + 1) + phiVal(i + 1, j) + phiVal(i,j) + phiVal(i, j + 1) 
		- (phiVal(i, j + 1) + phiVal(i,j) + phiVal(i - 1, j) + phiVal(i - 1, j + 1)));
	return Der;

};

double DERY_jMinusHalf(int i, int j){
	Der = (1 / ( 4 * dx_bar )) * 
		( phiVal(i + 1, j) + phiVal(i + 1, j - 1) + phiVal(i, j - 1) + phiVal(i, j) 
		- (phiVal(i, j) + phiVal(i,j - 1) + phiVal(i - 1, j - 1) + phiVal(i - 1, j)));
	return Der; 
};

double DERX(int i, int j){
	Der = (1/2) * (DERX_iPlusHalf + DERY_iMinusHalf);
	return Der;

};

double DERY(int i, int j){
	Der = (1/2) * (DERY_iPlusHalf + DERY_iMinusHalf);
	return Der;
}
/**************** MAG2 funcions ****************/
double MAG2_iPlusHalf(int i, int j){
	double mag, mag2;
	mag = pow(DERX_iPlusHalf(i, j), 2) 
			+ pow(DERY_iPlusHalf(i,j), 2);
	mag2 = pow(mag, 2);
};

double MAG2_iMinusHalf(int i, int j){
	double mag, mag2;
	mag = pow(DERX_iMinusHalf(i, j), 2) 
			+ pow(DERY_iMinusHalf(i,j), 2);
	mag2 = pow(mag, 2);
};
double MAG2_jPlusHalf(int i, int j){
	double mag, mag2;
	mag = pow(DERX_jPlusHalf(i, j), 2) 
			+ pow(DERY_jPlusHalf(i,j), 2);
	mag2 = pow(mag, 2);
};
double MAG2_iMinusHalf(int i, int j){
	double mag, mag2;
	mag = pow(DERX_iMinusHalf(i, j), 2) 
			+ pow(DERY_iMinusHalf(i,j), 2);
	mag2 = pow(mag, 2);
};


/**************** A funcions ****************/
double A_iPlusHalf(int i, int j){
	double A, num;
	num = pow(DERX_iPlusHalf(i,j) , 4) + pow(DERY_iPlusHalf(i,j), 4);
	A = a_s * (1 + epsilon_prime * ( num / MAG2_iPlusHalf(i,j) );
};

double A_iMinusHalf(int i, int j){
	double A, num;
	num = pow(DERX_iMinusHalf(i,j) , 4) + pow(DERY_iMinusHalf(i,j), 4);
	A = a_s * (1 + epsilon_prime * ( num / MAG2_iMinusHalf(i,j) );
};
double A_jPlusHalf(int i, int j){
	double A, num;
	num = pow(DERX_jPlusHalf(i,j) , 4) + pow(DERY_jPlusHalf(i,j), 4);
	A = a_s * (1 + epsilon_prime * ( num / MAG2_jPlusHalf(i,j) );
};

double A_jMinusHalf(int i, int j){
	double A, num;
	num = pow(DERX_jMinusHalf(i,j) , 4) + pow(DERY_jMinusHalf(i,j), 4);
	A = a_s * (1 + epsilon_prime * ( num / MAG2_jMinusHalf(i,j) );
};

double A(i,j){
	double A, num;
	num = pow(DERX(i,j) , 4) + pow(DERY(i,j), 4);
	A = a_s * (1 + epsilon_prime * ( num / MAG2(i,j) );
};

/**************** A' funcions ****************/
double Aprime_iPlusHalf(int i, int j){
	double A, num;
	num = pow(DERX_iPlusHalf(i,j), 2) - pow(DERY_iPlusHalf(i,j), 2);
	A = a12 * DERX_iPlusHalf(i,j) * DERY_iPlusHalf(i,j) * (num / MAG2_iMinusHalf(i,j));
	return A;
};

double Aprime_iMinusHalf(int i, int j){
	double A, num;
	num = pow(DERX_iMinusHalf(i,j), 2) - pow(DERY_iMinusHalf(i,j), 2);
	A = a12 * DERX_iMinusHalf(i,j) * DERY_iMinusHalf(i,j) * (num / MAG2_iMinusHalf(i,j));
	return A;
};
double Aprime_jPlusHalf(int i, int j){
	double A, num;
	num = pow(DERX_jPlusHalf(i,j), 2) - pow(DERY_jPlusHalf(i,j), 2);
	A = a12 * DERX_jPlusHalf(i,j) * DERY_jPlusHalf(i,j) * (num / MAG2_jMinusHalf(i,j));
	return A;
};

double Aprime_jMinusHalf(int i, int j){
	double A, num;
	num = pow(DERX_jMinusHalf(i,j), 2) - pow(DERY_jMinusHalf(i,j), 2);
	A = a12 * DERX_jMinusHalf(i,j) * DERY_jMinusHalf(i,j) * (num / MAG2_jMinusHalf(i,j));
	return A;
};

double A_prime(i,j){
	double A, num;
	num = pow(DERX(i,j), 2) - pow(DERY(i,j), 2);
	A = a12 * DERX(i,j) * DERY(i,j) * (num / MAG2(i,j));
	return A;
};

/*
"The arrays JR, JL, JT, JB respectively handle the gradient terms ("order parameter fluxes") from 
the phi equation on the right, left, top and bottom edges of the finite volume ventered around
the point (i,j).
*/
double JR(int i, int j){
	double jr;
	jr = A_iPlusHalf(i,j) * (A_iPlusHalf(i,j) * DERX_iPlusHalf(i,j)
				- Aprime_iPlusHalf(i,j) * DERY_iPlusHalf(i,j));
	return jr;
};

double JL(int i, int j){
	double jl;
	jl = A_iMinusHalf(i,j) * (A_iMinusHalf(i,j) * DERX_iMinusHalf(i,j)
				- Aprime_iMinusHalf(i,j) * DERY_iMinusHalf(i,j));
	return jl;
};

double JT(int i, int j){
	double jt;
	jt = A_jPlusHalf(i,j) * (A_jPlusHalf(i,j) * DERX_jPlusHalf(i,j)
				- Aprime_jPlusHalf(i,j) * DERY_jPlusHalf(i,j));
	return jt;

};

double JB(int i, int j){
	double jb;
	jb = A_jMinusHalf(i,j) * (A_jMinusHalf(i,j) * DERX_jMinusHalf(i,j)
				- Aprime_jMinusHalf(i,j) * DERY_jMinusHalf(i,j));
	return jb;
};

//phase field interpolation functions (also used in Karma and Rappel)
double g_Prime(int i, int j){
	double g;
	g = - phiVal(i,j) + pow( phiVal(i,j), 3);
	return g;
};

double P_prime(int i, int j){
	double p, p2;
	p = 1 - pow(phiVal(i,j), 2);
	p2 = pow(phiVal(i,j), 2);
	return p2;
};

double h_Prime(i,j){
	reutrn (1/2);
};

double UpdatePHI(int i, int j){
	double newPHI;
	newPHI = phiVal(i,j) + (dt_bar / pow)

};

//evaluates the laplacian for phi(i,j). Eqn B.5 in Nik's book may be used as a more isotropic laplacian
//for PERIODIC BOUNDARY CONDITIONS - right now doesn't deal with 0 nodes well...
double Laplacian_phi(int i, int j){

	double laplacian_ij; 

	laplacian_ij = phiVal(i + 1 , j) + phiVal(i - 1, j)
			+ phiVal( i, j + 1 ) + phiVal(i, j - 1) 
			- 4 * ( phiVal(i,j) );
	
	return laplacian_ij;
};

double MU(int i, int j){
	double mu = - ( W / pow(dx, 2) ) * Laplacian_phi(i,j) + df(i,j);
	return mu;
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

double timeMarchC(int i, int j){

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
	double Temp, Tm, L, c_P, alpha, tao;

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
	W = 0.25;
	M = 1; //lookup realistic value for M.

	//model C specific variables
	epsilon_4 = 1; //anisotropy parameter form equation 5.25
	epsilon_prime = 4 * epsilon_4 / a_s ;
	a_s = 1 - 3 * epsilon_4;
	a12 = 4 * a_s * epsilon_prime;

	Temp = 20;
	Tm = 20;
	c_P = 1;
	L = 22;
	alpha = 5;
	tao = 3;
	lambda_bar = 4;

	U = (Temp - Tm)* c_P / (L);
	D_bar = (alpha * tao) / pow(W, 2);

	dt_bar = dt / tao;
	dx_bar = dx / pow(0.25, 1/2);


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
				phi[i][j] = timeMarchC(i,j);
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