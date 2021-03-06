// compilacion: gcc rastreging.c -o out -lm 
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#define POP_SIZE 2048
#define IND_SIZE 64
#define STEP 1024 // 1/STEP double precission
// bits = log2 STEP 
//16 - bits = numBits domain size representation
#define BITS 10 
#define ITERATIONS 40000

// domain intervale limit = pow (2, ind_size/(num_vars) - bits ) /2  

double Hist_res[ITERATIONS]; // historico de selecciones fitness

struct Individuo
{
	char genes[IND_SIZE];
};

// individuo poblacion inicial

struct Individuo pop[POP_SIZE]; 
int fittest = 0; 
int secondFittest = 0; 
int mutate_var = 0; // no booleans on C :0
double res = 0.0;




// tamaño de representacion string binaria depende del tamaño de busqueda del dominio del problema
// 2^|log2(size(dominio))|
// -1000 a 1000
// 0000 = -1000
// 1111 = 1000
// size = 2000
// 2^11 = 2048

void generatePopulation();
void evaluate();
void cruce(); 
void mutate(); 
void imprimir_poblacion();
void imprimir_individuo();
void decodificar_individuo();

//https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3



int main()
{

time_t t;
srand((unsigned int) time(&t));

generatePopulation();
printf("Población Inicial:\n");
for (int i = 0; i<ITERATIONS; i++)
{
imprimir_poblacion();
evaluate();
Hist_res[i] = res;
cruce();
if (mutate_var)
{
//	printf("var %d en interacion %d\n", mutate_var, i);
	mutate();
		

} 
//printf("\nLa poblacion %d es\n", i+1);
//imprimir_poblacion();
//evaluate(); // final sol
}
for (int i = 0; i<ITERATIONS; i++)
{
	printf("it %d\t%f\n", i, Hist_res[i]);
}

imprimir_poblacion();
evaluate(); // final sol
//printf("%f", Hist_res[ITERATIONS-1]);


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
/*			printf("%d", s);
			printf("xx%d", pop[j].genes[i]);*/
		}
	//	printf("\n");
	}

}


//res+= (p[i]*p[i]) - (10 * cos(p[i])) + 10; rastreging
// f(x,y) = x^2 = 2y^2 - 0.3cos(3pix) - 0.4cos(4piy) + 0.7
// f(x,y,w,z) = 100(y - x^2)^2 + (1-x)^2 +90(z-w^2)^2 + (1 - w)^2 +10.1((y-1)^2 + (z-1)^2)

void evaluate()
{

	 
	float aux_fittest[POP_SIZE];
	
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
	
	//asigna nueva poblacion en estructura de datos
	memcpy(&pop, &new_pop, sizeof new_pop);
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
{
	for ( int i = 1; i < POP_SIZE; i++)
	{//numvars
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
