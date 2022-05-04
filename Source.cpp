/*Samuel Celli
* qr3428
* Genetic algorithm coding assignment
* CS461 Moayed Daneshyari
*/


#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <time.h>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

const int CHROM_SIZE = 16;
const int POP_SIZE = 8;
const int GEN_NUM = 200;
const double CROSSOVER_PROB = 0.7;
const double MUTATION_PROB = 0.01;

//The function to find a maximum for
double func(double x, double y) {
	double total = (pow((1.0 - x), 2.0) * pow(M_E, (-pow(x, 2) - pow(y + 1, 2))))
		- ((x - pow(x, 3) - pow(y, 3)) * pow(M_E, (-pow(x, 2) - pow(y, 2))));
	return total;
}

//conversion of int value of genes into function domain
double real_input_val(double x) {
	return (x * double(4.0 / 255.0) - 2);
}

/*This class holds the bit arrays for the chromosones and handles
* all operations on them
*/
class Chromosone {
public:
	Chromosone();
	Chromosone(int);
	Chromosone(bool*);
	int getX();
	int getY();
	static int get_mut_count() { return mutation_count; }
	bool* get_genes();
	double get_fitness(double (double, double));
	static bool * swap(Chromosone*, Chromosone*, int);
	void mutate(int);
private:
	bool genes[CHROM_SIZE];
	static int mutation_count;
};

//For debugging and print function
int Chromosone::mutation_count = 0;

Chromosone::Chromosone() {
	int init = rand() % (int)(pow( 2.0, (double)CHROM_SIZE)-1);
	for (int i = 0; i < CHROM_SIZE; i++) {
		this->genes[i] = init % 2;
		init = init / 2;
	}
	return;
}

Chromosone::Chromosone(int init) {
	for (int i = 0; i < CHROM_SIZE; i++) {
		this->genes[i] = init % 2;
		init = init / 2;
	}
	return;
}

Chromosone::Chromosone(bool *init_genes) {
	for (int i = 0; i < CHROM_SIZE; i++) {
		this->genes[i] = init_genes[i];
	}
	return;
}

//converts the first half of bits into an int
int Chromosone::getX() {
	int total = 0;
	for (int i = 0; i < CHROM_SIZE/2; i++) {
		total *= 2;
		total += (int)genes[i];
	}
	return total;
}

//converts the second half of bits into an int
int Chromosone::getY() {
	int total = 0;
	for (int i = CHROM_SIZE/2; i < CHROM_SIZE; i++) {
		total *= 2;
		total += (int)genes[i];
	}
	return total;
}

//returns a copy of the array of bits 
bool* Chromosone::get_genes() {
	static bool temp[CHROM_SIZE] = {};
	for (int i = 0; i < CHROM_SIZE; i++) {
		temp[i] = this->genes[i];
	}
	return temp;
}

/* Evaluates the chromosones against the function to determine fitness
* since we are looking for a maximum and the function does return positive
* values, negative values are discarded with a fitness of 0
*/
double Chromosone::get_fitness(double foo(double, double)) {
	double x = real_input_val((double)this->getX());
	double y = real_input_val((double)this->getY());
	double result = foo(x, y);
	if (result < 0)
		result = 0;
	return result;
}

/*returns a set of chromosones based on 2 parents and a breakpoint
*/
bool * Chromosone::swap(Chromosone *Chrom1, Chromosone *Chrom2, int pos) {
	static bool temp[CHROM_SIZE];
	for (int i = 0; i < pos; i++) {
		temp[i] = Chrom1->genes[i];
	}
	for (int i = pos; i < CHROM_SIZE; i++) {
		temp[i] = Chrom2->genes[i];
	}
	return temp;
}

/*mutates a bit at the position provided or a random one if not provided
*/
void Chromosone::mutate(int pos = CHROM_SIZE) {
	if (pos >= CHROM_SIZE || pos < 0) {
		pos = rand() % CHROM_SIZE;
	}
	this->genes[pos] = !this->genes[pos];
	mutation_count++;
	return;
}

/*This class contains the population of chromosones and all operations with them
* including handling fitness and fitness ratio, as well as creating children from
* parents. It also prints data and keeps track of the amount of generations.
*/
class Generation {
public:
	Generation();
	~Generation();
	int get_gen() { return current_gen; }
	int get_pop() { return pop; }
	void print_fitness();
	Chromosone* get_parent();
	void calc_pop_fitness();
	void calc_fitness_ratio();
	void random_gen();
	void set_child(Chromosone**);
	void det_mutation(Chromosone*);
	Chromosone ** create_child();

private:
	int pop;
	static int current_gen;
	static int crossover_count;
	Chromosone* population[POP_SIZE] = {};
	double pop_fitness[POP_SIZE] = {};
	double fitness_ratio[POP_SIZE] = {};
};

//Used for debugging and printing output
int Generation::current_gen = 0;
int Generation::crossover_count = 0;

Generation::Generation() {
	pop = 0;
	current_gen++;
	return;
}

Generation::~Generation() {
	for (int i = 0; i < pop; i++)
		delete population[i];
	return;
}

//prints a report with fitness and input values for a generation
void Generation::print_fitness() {
	cout << "******  Fitness Report  ********" << endl;
	for (int i = 0; i < pop; i++) {
		cout << "POP# " << setw(10) << std::left << i 
			<< setw(10) << setprecision(5) << pop_fitness[i] 
			<< setw(10) << setprecision(5) << real_input_val(population[i]->getX())
			<< setw(10) << setprecision(5) << real_input_val(population[i]->getY())
			<< endl;
	}
	cout << "There were " << current_gen << " generations" << endl;
	cout << "There were " << Chromosone::get_mut_count() << " mutations" << endl;
	cout << "There were " << crossover_count <<" crossovers" << endl;
}

/*Chooses a Chromosone based on its relative fitness and returns a pointer to it
*/
Chromosone* Generation::get_parent() {
	double val = (double)rand() / RAND_MAX;
	Chromosone* parent=NULL;
	for (int y = 0; y < pop; y++) {
		if (val < fitness_ratio[y]) {
			parent = population[y];
			break;
		}
	}
	if (!parent)
		parent = population[POP_SIZE - 1];
	return parent;
}

/*Populates the array of fitness values
*/
void Generation::calc_pop_fitness() {
	for (int i = 0; i < pop; i++) {
		pop_fitness[i] = population[i]->get_fitness(func);
	}
	return;
}

/*Populates the array of fitness ratio values for selecting random chromosones
* to have children
*/
void Generation::calc_fitness_ratio() {
	double total = 0;
	for (int i = 0; i < pop; i++) {
		total += pop_fitness[i];
	}
	for (int i = 0; i < pop; i++) {
		fitness_ratio[i] += pop_fitness[i] / total;
		if (i + 1 < pop) {
			fitness_ratio[i + 1] = fitness_ratio[i];
		}
	}
	return;
}

/*Generates a completely random generation to start
*/
void Generation::random_gen() {
	for (pop; pop < POP_SIZE; pop++) {
		population[pop] = new Chromosone;
	}
	calc_pop_fitness();
	calc_fitness_ratio();
}

/*places newly created children into a new generation
*/
void Generation::set_child(Chromosone ** new_children) {
	if (pop >= POP_SIZE - 1)
		return;
	population[pop] = new_children[0];
	pop++;
	population[pop] = new_children[1];
	pop++;
}

/*determines if a chromosone should have a mutation
*/
void Generation::det_mutation(Chromosone* chrom) {
	for (int i = 0; i < CHROM_SIZE; i++) {
		if (rand() % 1000 < MUTATION_PROB * 1000)
			chrom->mutate(i);
	}
	return;
}

/* Gets 2 unique parents, decides if they should crossover, generates children,
* subjects children to mutation, and returns an array of 2 new Chromosones
*/
Chromosone** Generation::create_child() {
	Chromosone* parents[2] = {};
	static Chromosone * children[2];
	parents[0] = get_parent();
	parents[1] = get_parent();

	while (parents[0] == parents[1])
		parents[1] = get_parent();

	if (rand() % 100 < 100 * CROSSOVER_PROB) {
		int pos = 1 + rand() % CHROM_SIZE-2;
		children[0] = new Chromosone(Chromosone::swap(parents[0], parents[1], pos));
		children[1] = new Chromosone(Chromosone::swap(parents[1], parents[0], pos));
		crossover_count++;
	}
	else {
		children[0] = new Chromosone(parents[0]->get_genes());
		children[1] = new Chromosone(parents[1]->get_genes());
	}
	det_mutation(children[0]);
	det_mutation(children[1]);
	return children;
}

int main() {
	srand(time(NULL));
	Generation* Current;
	Generation* Newest;
	Current = new Generation;
	Current->random_gen();
	Current->print_fitness();
	while (Current->get_gen() < GEN_NUM) {
		Newest = new Generation;
		for (int i=0; i<POP_SIZE; i+=2)
			Newest->set_child(Current->create_child());
		Newest->calc_pop_fitness();
		Newest->calc_fitness_ratio();
		delete Current;
		Current = Newest;
		Newest = NULL;
	}
	Current->print_fitness();
}