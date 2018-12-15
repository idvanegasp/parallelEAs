// compilacion: gcc rastreging.c -o out -lm -lpthread
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
#define NUM_THREADS 4

// numero de individuos a evaluar por hilo POP_SIZE/NUM_THREADS

double Hist_res[ITERATIONS]; // historico de selecciones fitness

struct Individuo
{
	char genes[IND_SIZE];
};


// individuo poblacion inicial
struct Individuo pop[POP_SIZE]; 

// variables de estado, metricas
int fittest = 0; 
int secondFittest = 0; 
double aux_fittest[POP_SIZE];
int mutate_var = 0; // no booleans on C
double res = 0.0;
int iterator = 0;
void* generatePopulation(void* arg);
void* evaluate(void* arg); 
void do_selection();
void cruce(); 
void mutate(); 
void imprimir_poblacion();
void imprimir_individuo();
void decodificar_individuo();

//https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3


int main()
{

pthread_t threads[NUM_THREADS];

time_t t;
srand((unsigned int) time(&t));


// generacion paralela de poblacion inicial
for (int l = 0; l<NUM_THREADS; l++)
{
double *p;
pthread_create(&threads[l], NULL, generatePopulation, (void*)(p));
}

for (int l = 0; l<NUM_THREADS; l++)
{
pthread_join(threads[l], NULL);
}
printf("Población Inicial:\n");

// flujo de evolución
for (int i = 0; i<ITERATIONS; i++)
{
imprimir_poblacion();

// Evaluacion paralela de generaciones
iterator = 0;
for (int l = 0; l<NUM_THREADS; l++)
{
double *p;
pthread_create(&threads[l], NULL, evaluate, (void*)(p));
}

for (int l = 0; l<NUM_THREADS; l++)
{
pthread_join(threads[l], NULL);
}

/*
for (int l = 0; l<POP_SIZE; l++)
{
	printf("\nit %d\t%f", l, aux_fittest[l]);
}
*/

do_selection();
//actualiza historico para llevar metricas de fitness
// (fuera de alcance)
Hist_res[i] = res;
cruce();
// evalua probabilidad de mutacion
if (mutate_var)	mutate();

//printf("Generacion %d:\n", i);
}


imprimir_poblacion();
double *p;
evaluate(p); // final sol
//printf("%f", Hist_res[ITERATIONS-1]);
do_selection();

//imprime evolucion del fitness de cada generacion
/*for (int i = 0; i<ITERATIONS; i++)
	printf("\nit %d\t%f", i, Hist_res[i]);
*/
return 0;
}




/*
Function SET evaluation chooser
Optimization type min max
int fn;
scanf("%d", &fn);


float *first ;

first = malloc( 2 * sizeof(float));
first[0] = (float)(rand() & 0xFF) / 10.0f;
first[1] = (float)(rand() & 0xFF) / 10.0f; 
*/




void* generatePopulation(void* arg){

	int core_thread = iterator++;

	for (int i = core_thread*POP_SIZE/NUM_THREADS; i<(core_thread+1)*POP_SIZE/NUM_THREADS; i++)
	{
		printf("\nHilo %d, indice %d", core_thread, i);
		strcpy( pop[i].genes, "\0");
		int s; // random bit
		for (int j = 0; j < IND_SIZE; j++)
		{
			s = rand()%2;
			pop[i].genes[j] = (char)s;
		}
	}
}


//res+= (p[i]*p[i]) - (10 * cos(p[i])) + 10; rastreging
// f(x,y) = x^2 = 2y^2 - 0.3cos(3pix) - 0.4cos(4piy) + 0.7
// f(x,y,w,z) = 100(y - x^2)^2 + (1-x)^2 +90(z-w^2)^2 + (1 - w)^2 +10.1((y-1)^2 + (z-1)^2)

// dom intervale limit = pow (2, ind_size/(num_vars) - bits ) /2  
void* evaluate(void* arg)
{
	int core_thread = iterator++;

	
	for (int j = core_thread*POP_SIZE/NUM_THREADS; j<(core_thread+1)*POP_SIZE/NUM_THREADS; j++)
	{	
		double fit_val = 0.0; 
		double aux_X = 0.0;
		double aux_Y = 0.0;
		double aux_W = 0.0;
		double aux_Z = 0.0;

/// mejorar: X_Y_W_Z[] = funcion_decodificarIndividuo(chrom, res[])
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
		

		fit_val = 100*pow((aux_Y - pow(aux_X,2)), 2) + pow((1-aux_X),2) + 90*pow((aux_Z-pow(aux_W,2)),2) + pow((1 - aux_W),2) +10.1*(pow((aux_Y-1),2) + pow((aux_Z-1),2));
		//printf("X%d : %f, Y%d : %f;  W%d : %f;  Z%d : %f; f%d = %f\n", j+1,aux_X, j+1, aux_Y, j+1, aux_W, j+1, aux_Z, j+1, fit_val);
		aux_fittest[j] = fit_val; 
	}	
}


void do_selection()
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
	printf("El individuo minimo es %d y el segundo %d\n", fittest+1, secondFittest+1);
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
	
/*
	for (int j = 0; j < POP_SIZE; j++)
	{
		for (int i = 0; i < IND_SIZE; i++)
		{
			printf("%d", new_pop[j].genes[i]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
*/	memcpy(&pop, &new_pop, sizeof new_pop);
}


void no_mutate()
{
	for ( int i = 0; i < POP_SIZE; i++)
	{//numvars
	int mut = rand()%2;
	if (mut)
	{
		int k = rand()%IND_SIZE ;
		if(pop[i].genes[k]%2==1)
		{
			pop[i].genes[k] = 0;
	//			printf("Muta Ind%d cromosoma%d\n", i+1, k+1);
		}else pop[i].genes[k] = 1;
	}
	
	// exchange rows
	}
	mutate_var = 0;
/*	for ( int i = IND_SIZE/2; i <= (IND_SIZE/2)+k; i++)
	{

	}
*/
}

void mutate()
{	// indice 1 evita perder el mejor individuo de todas las generaciones
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
/*	for ( int i = IND_SIZE/2; i <= (IND_SIZE/2)+k; i++)
	{

	}
*/
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
//	printf("\n");
//	printf("\n");
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

// coding challenge 35 genetic algorithm
