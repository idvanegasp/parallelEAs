#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <iostream>

//using namespace std;

#define POP_SIZE 4
#define IND_SIZE 32
#define STEP 64 
// bits = log STEP
// dom intervale limit = pow (2, ind_size/2 - bits ) /2  
#define ITERATIONS 1000

float res = 0;


// compilacion: gcc rastreging.c -o out -lm ; gracias Cynthia

// tamaño de representacion string binaria depende del tamaño de busqueda del dominio del problema
// 2^|log2(size(dominio))|
// -1000 a 1000
// 0000 = -1000
// 1111 = 1000
// size = 2000
// 2^11 = 2048


// size * 2 x,y,
// valor real = bitstring - size/2

void generatePopulation();
void evaluate();
void cruiceBestPop(); 
void mutateInd(); // 1 in a Thousand or more; low prob. 
void imprimir_poblacion();

//https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3


struct Individuo
{
	//size_t ind_size;
	char genes[IND_SIZE];
};


// struct ind gen1 gen2 ... gen11

// individuo poblacion inicial



struct Individuo pop[POP_SIZE]; 
int fittest = 0; 
int secondFittest = 0; 


int main()
{

time_t t;
srand((unsigned int) time(&t));

generatePopulation();

imprimir_poblacion();

float *first ;

first = malloc( 2 * sizeof(float));
first[0] = (float)(rand() & 0xFF) / 10.0f;
first[1] = (float)(rand() & 0xFF) / 10.0f; 

/*
Function SET evaluation chooser

int fn;
scanf("%d", &fn);
*/


evaluate();

imprimir_poblacion();

return 0;
}


void generatePopulation(){

	strcpy( pop[0].genes, "\0");
	strcpy( pop[1].genes, "\0");
	strcpy( pop[2].genes, "\0");
	strcpy( pop[3].genes, "\0");
	int s;
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


//res+= (p[i]*p[i]) - (10 * cos(p[i])) + 10;
void evaluate()
{
	double fit_val = 0.0; 
	for (int j = 0; j < POP_SIZE; j++)
	{
		double aux_X = 0.0;
		double aux_Y = 0.0;
		// 6 = log2 (step = 64)
		for (int i = 0; i < IND_SIZE/2; i++)
		{
			aux_X += pop[j].genes[i]*pow(2 , ((IND_SIZE/2)-i)-1); 

		}
		aux_X = (aux_X/STEP) - pow(2,(IND_SIZE/2)-6)/2;
		for (int i = IND_SIZE/2; i < IND_SIZE; i++)
		{
			aux_Y += pop[j].genes[i]*pow(2 , (IND_SIZE-i)-1); 

		}
		aux_Y = (aux_Y/STEP) - pow(2,(IND_SIZE/2)-6)/2;
		printf("Ind X%d, vale %f\n", j+1, aux_X);
		printf("Ind Y%d, vale %f\n", j+1, aux_Y);
	}
}


void mutate()
{
	int k = rand()% (IND_SIZE/2) ;
	for ( int i = 0; i <= k; i++)
	{
		// exchange rows
	}
	for ( int i = IND_SIZE/2; i <= (IND_SIZE/2)+k; i++)
	{

	}
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
