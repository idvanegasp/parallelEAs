#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <iostream>

//using namespace std;

#define POP_SIZE 4
#define ITERATIONS 100

float res = 0;


// compilacion: gcc rastreging.c -o out -lm ; gracias Cynthia

// tamaño de representacion string binaria depende del tamaño de busqueda del dominio del problema
// 2^|log2(size(dominio))|
// -1000 a 1000
// 0000 = -1000
// 1111 = 1000
// size = 2000
// 2^11 = 2048

void evaluate(float p[]);
void generatePopulation();
void cruiceBestPop(float p[]); 
void mutateInd(); // 1 in a Thousand or more; low prob. 
void imprimir_poblacion();

//https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3



struct Individuo
{
	//size_t ind_size;
	char genes[11];
};


// struct ind gen1 gen2 ... gen11

// individuo poblacion inicial

size_t ind_size = 12;

struct Individuo pop[POP_SIZE]; 



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

evaluate(first);

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
		for (int i = 0; i < ind_size; i++)
		{
			s = rand()%2;
			pop[j].genes[i] = (char)s;
/*			printf("%d", s);
			printf("xx%d", pop[j].genes[i]);*/
		}
	//	printf("\n");
	}

}

void evaluate(float p[])
{
for (int i = 0; i<2; i++)
	res+= (p[i]*p[i]) - (10 * cos(p[i])) + 10;
{
}
printf("%.2f %.2f\n", p[0], p[1]);
printf(" f: %.2f\n", res);
}

void imprimir_poblacion()
{
	for (int j = 0; j < POP_SIZE; j++)
	{
		for (int i = 0; i < ind_size; i++)
		{
			printf("%d", pop[j].genes[i]);
		}
		printf("\n");
	}

	printf("\n");
	printf("\n");
}
