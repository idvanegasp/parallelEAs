#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define POP_SIZE 4
#define ITERATIONS 100

float res = 0;

// compilacion: gcc rastreging.c -o out -lm ; gracias Cynthia

// tamaño de representacion string binaria depende del tamaño de busqueda del dominio del problema
// 2^|log2(size(dominio))|

void evaluate(float p[]);
void generatePopulation();
void cruiceBestPop(float p[]); 
void mutateInd(); // 1 in a Thousand or more; low prob. 

int main()
{
// individuo poblacion inicial
time_t t;
srand((unsigned int) time(&t));
float *first ;

first = malloc( 2 * sizeof(float));
first[0] = (float)(rand() & 0xFF) / 10.0f;
first[1] = (float)(rand() & 0xFF) / 10.0f; 

evaluate(first);

return 0;
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

void generatePop();
