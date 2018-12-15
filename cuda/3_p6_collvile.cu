// compilacion: nvcc 3_p6_collvile.cu -o collville

#include <cuda_runtime.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#define POP_SIZE 2048
#define IND_SIZE 64
#define STEP 1024 
// bits = log2 STEP 
//16 - bits = numBits domain size representation
#define BITS 10
#define ITERATIONS 40000

double Hist_res[ITERATIONS];

struct Individuo
{
	char genes[IND_SIZE];
};


struct Individuo pop[POP_SIZE]; 

char* poblacion;
int fittest = 0; 
int secondFittest = 0; 
double aux_fittest[POP_SIZE];
int mutate_var = 0; // no booleans on C :0
double res = 0.0;


void generatePopulation();
void evaluate();
void select();
void cruce(); 
void mutate(); 
void imprimir_poblacion();
void imprimir_individuo(int idx);
void decodificar_individuo(int idx);

__global__ void evaluateKernel( char* pop, double* aux_fittest_gpu, int chrom_size, int pop_size, int step, int bits ){

	int h = threadIdx.x + blockDim.x*blockIdx.x ;
	///int h2 = threadIdx.y+ blockDim.y*blockIdx.y ;
	double fit_val = 0.0;
	double aux_X = 0.0;
	double aux_Y = 0.0;
	double aux_W = 0.0;
	double aux_Z = 0.0;

if (h < pop_size)
{
	for (int k = 0; k<chrom_size/4; k++)
		aux_X = aux_X+ pop[h]*powf(2 , ((chrom_size/4)-k)-1);
	aux_X < (aux_X/step) - powf(2,(chrom_size/4)-bits)/2;

	for (int k = chrom_size/4; k<chrom_size/2; k++)
		aux_Y = aux_Y+ pop[h]*powf(2 , ((chrom_size/2)-k)-1);
		aux_Y = (aux_Y/step) - powf(2,(chrom_size/4)-bits)/2;
	for (int k = chrom_size/2; k<chrom_size - 16; k++)
		aux_W = aux_W +pop[h]*powf(2 , ((chrom_size-16)-k)-1);
		aux_W = (aux_W/step) - powf(2,(chrom_size/4)-bits)/2;
	for (int k = (chrom_size-16); k<chrom_size; k++)
		aux_Z = aux_Z+ pop[h]*powf(2 , (chrom_size-k)-1);
		aux_Z = (aux_Z/step) - powf(2,(chrom_size/4)-bits)/2;

		fit_val = 100*powf((aux_Y - powf(aux_X,2)), 2) + powf((1-aux_X),2) + 90*powf((aux_Z-powf(aux_W,2)),2) + powf((1 - aux_W),2) +10.1*(powf((aux_Y-1),2) + powf((aux_Z-1),2));

		aux_fittest_gpu[h] = fit_val; 
}
	
}



//https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3

// Inicio bloque principal

int main(void )
{

time_t t;
srand((unsigned int) time(&t));
cudaDeviceProp deviceProp;
cudaGetDeviceProperties(&deviceProp, 0);
generatePopulation();
printf("Población Inicial:\n");
//imprimir_poblacion();
printf("starting evolution ... %s" , deviceProp.name);

	cudaSetDevice(0);
	// allocate device memory
	// allocat device memory
	size_t nBytes = sizeof(char *)*POP_SIZE*IND_SIZE;
	char* gpuRef; 
	gpuRef = (char*)malloc(nBytes);
	char* d_P;
	double* aux_fittest_gpu;// = (double *)malloc(POP_SIZE);
	poblacion = (char*)malloc(nBytes);
	memset(gpuRef,0, nBytes);
//	memset(aux_fittest_gpu, 0, POP_SIZE);
	cudaMalloc((char** ) &d_P, nBytes);
	cudaMalloc((double**) &aux_fittest_gpu, POP_SIZE);
	

	dim3 blockDim(32,32);
	dim3 gridDim(POP_SIZE/32,IND_SIZE/32);


for (int i = 0; i<ITERATIONS; i++)
{
cudaMemcpy(d_P, poblacion, nBytes, cudaMemcpyHostToDevice);
evaluateKernel<<<gridDim, blockDim>>>(d_P, aux_fittest_gpu, IND_SIZE, POP_SIZE, STEP, BITS);
cudaDeviceSynchronize();
cudaMemcpy(gpuRef, d_P, nBytes, cudaMemcpyDeviceToHost);
cudaMemcpy(aux_fittest, aux_fittest_gpu, nBytes, cudaMemcpyDeviceToHost);
select();
Hist_res[i] = res;
cruce();
//imprimir_poblacion();
if (mutate_var)	mutate();
		
}
evaluate();
select();
//imprime evolucion del fitness de cada generacion
/*for (int i = 0; i<ITERATIONS; i++)
	printf("\nit %d\t%f", i, Hist_res[i]);
*/
cudaFree(d_P);
cudaFree(aux_fittest_gpu);

return 0;
}

// Fin bloque principal

void generatePopulation(){

	for (int i = 0; i<POP_SIZE; i++)
	{
		strcpy( pop[i].genes, "\0");
	}
	int s; // random bit
	for (int j = 0; j < POP_SIZE; j++)
	{	
		for (int i = 0; i < IND_SIZE; i++)
		{
			s = rand()%2;
			pop[j].genes[i] = (char)s;
		}
	}

}


//res+= (p[i]*p[i]) - (10 * cos(p[i])) + 10; rastreging
// f(x,y) = x^2 = 2y^2 - 0.3cos(3pix) - 0.4cos(4piy) + 0.7
// f(x,y,w,z) = 100(y - x^2)^2 + (1-x)^2 +90(z-w^2)^2 + (1 - w)^2 +10.1((y-1)^2 + (z-1)^2)

void evaluate()
{

for (int j = 0; j < POP_SIZE; j++)
	{	
		double fit_val = 0.0;
		double aux_X = 0.0;
		double aux_Y = 0.0;
		double aux_W = 0.0;
		double aux_Z = 0.0;
		// bloques altamente paralelizables
		for (int i = 0; i < IND_SIZE/4; i++)
		{
			aux_X += pop[j].genes[i]*pow(2 , ((IND_SIZE/4)-i)-1); 

		}
		aux_X = (aux_X/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
		for (int i = IND_SIZE/4; i < IND_SIZE/2; i++)
		{
			aux_Y += pop[j].genes[i]*pow(2 , ((IND_SIZE/2)-i)-1); 

		}
		aux_Y = (aux_Y/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
		for (int i = IND_SIZE/2; i < IND_SIZE-16; i++)
		{
			aux_W += pop[j].genes[i]*pow(2 , (IND_SIZE-16-i)-1); 

		}
		aux_W = (aux_W/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
		for (int i = IND_SIZE-16 ; i < IND_SIZE; i++)
		{
			aux_Z += pop[j].genes[i]*pow(2 , (IND_SIZE-i)-1); 

		}
		aux_Z = (aux_Z/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
		
		//res += (pow(aux_X, 2)) - (10 * cos(2*aux_X)) + 10;
		//res += (pow(aux_Y, 2)) - (10 * cos(2*aux_Y)) + 10;
		
//		res = pow(x, 2) + 2*pow(y, 2) -0.3*cos(3*pi*x) -0.4*cos(4*pi*y) + 0.7;
		fit_val = 100*pow((aux_Y - pow(aux_X,2)), 2) + pow((1-aux_X),2) + 90*pow((aux_Z-pow(aux_W,2)),2) + pow((1 - aux_W),2) +10.1*(pow((aux_Y-1),2) + pow((aux_Z-1),2));
		//printf("X%d : %f, Y%d : %f;  W%d : %f;  Z%d : %f; f%d = %f\n", j+1,aux_X, j+1, aux_Y, j+1, aux_W, j+1, aux_Z, j+1, res);
		aux_fittest[j] = fit_val; 
	}


	
}


void select()
{
	res = 0.0;
	int aux_fit_1 = 0;
	for (int j = 1; j < POP_SIZE; j++)
	{
		if (aux_fittest[j] < aux_fittest[aux_fit_1]){
			aux_fit_1 = j;		
		}
	}
	fittest = aux_fit_1;
	int aux_fit_2 = 0;
	if (aux_fit_1 == 0) aux_fit_2++;
	for (int j = aux_fit_2+1; j < POP_SIZE; j++)
	{
		if (j != aux_fit_1 && aux_fittest[j] < aux_fittest[aux_fit_2]){
			aux_fit_2 = j;		
		}
	}	
	secondFittest = aux_fit_2;
	if(abs(aux_fittest[aux_fit_1]-aux_fittest[aux_fit_2]) < 0.00001) mutate_var = 1;
	//printf("El individuo minimo es %d y el segundo %d\n", fittest+1, secondFittest+1);
	res=aux_fittest[fittest];
	decodificar_individuo(fittest);
}
void cruce()
{

	int s;
	s = rand()%IND_SIZE;
	
	char aux_ind_1[IND_SIZE];
	char aux_ind_2[IND_SIZE];

//	strcpy( aux_ind_1, &pop[fittest].genes);
//	strcpy( aux_ind_2, &pop[secondFittest].genes);
	
	for (int i = 0; i < IND_SIZE; i++)
	{
		aux_ind_1[i] = pop[fittest].genes[i];
		aux_ind_2[i] = pop[secondFittest].genes[i];
	}
	printf("\n");

	// crea los dos nuevos individuos del offspring y actualiza la poblacion
	struct Individuo new_pop[POP_SIZE]; 
	for (int i = 0; i<POP_SIZE; i++)
	{
		strcpy( new_pop[i].genes, "\0");
	}

// crea nueva poblacion
	for (int j = 0; j<POP_SIZE; j+=2)
	{
		for (int i = 0; i < IND_SIZE; i++)
		{
			new_pop[j].genes[i] = aux_ind_1[i];
			new_pop[j+1].genes[i] = aux_ind_2[i];
		}
	}
// actualiza la info de los hijos
	for (int j = 2; j<POP_SIZE; j+=2)
	{
		s = rand()%IND_SIZE;
		for (int i = s+1; i < IND_SIZE; i++)
		{
			new_pop[j].genes[i] = aux_ind_2[i];
			new_pop[j+1].genes[i] = aux_ind_1[i];
		}
	}
	memcpy(&pop, &new_pop, sizeof new_pop);
}


void mutate()
{
	for ( int i = 1; i < POP_SIZE; i++)
	{

	//numvars
	for ( int j = 0 ; j< IND_SIZE/16 ; j++)

	{
		int mut = rand()%2;
		if (mut)
		{
			int k = rand()%IND_SIZE/4 ;
			if(pop[i].genes[(j*IND_SIZE/4) + k]%2==1)
			{
				pop[i].genes[(j*IND_SIZE/4) + k] = 0;
	//			printf("Muta Ind%d cromosoma%d\n", i+1, k+1);
			}else pop[i].genes[(j*IND_SIZE/4) + k] = 1;
		}
	
	}

		// exchange rows
	}
	mutate_var = 0;
}

void imprimir_individuo(int k)
{
	for (int i = 0; i < IND_SIZE; i++)
	{
		printf("%d", pop[k].genes[i]);
	}
	printf("\n");
}
void imprimir_poblacion()
{
	for (int j = 0; j < POP_SIZE; j++)
	{
		for (int i = 0; i < IND_SIZE; i++)
		{
			printf("%d", pop[j].genes[i]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
}

void decodificar_individuo(int k)
{
	imprimir_individuo(k);
	double aux_X = 0.0;
	double aux_Y = 0.0;
	double aux_W = 0.0;
	double aux_Z = 0.0;

/// mejorar: X_Y_W_Z[] = funcion_decodificarIndividuo(chrom, res[])
	for (int i = 0; i < IND_SIZE/4; i++)
	{
		aux_X += pop[k].genes[i]*pow(2 , ((IND_SIZE/4)-i)-1); 
	}
	aux_X = (aux_X/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
	for (int i = IND_SIZE/4; i < IND_SIZE/2; i++)
	{
		aux_Y += pop[k].genes[i]*pow(2 , ((IND_SIZE/2)-i)-1); 
	}
	aux_Y = (aux_Y/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
	for (int i = IND_SIZE/2; i < IND_SIZE-16; i++)
	{
		aux_W += pop[k].genes[i]*pow(2 , (IND_SIZE-16-i)-1); 
		}
	aux_W = (aux_W/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
	for (int i = IND_SIZE-16 ; i < IND_SIZE; i++)
	{
		aux_Z += pop[k].genes[i]*pow(2 , (IND_SIZE-i)-1); 
	}
	aux_Z = (aux_Z/STEP) - pow(2,(IND_SIZE/4)-BITS)/2;
	printf("X: %f Y: %f W: %f Z: %f\n fitness %f\n", aux_X, aux_Y, aux_W, aux_Z, res);
}
