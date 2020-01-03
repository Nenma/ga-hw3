#include <iostream>
#include <math.h> 
#include <random>
#include <chrono> 
#include <iomanip> 
#define MIN -9223372036854775
#define MAX 9223372036854775
#define BITMAX 1000000
#define POPSIZE 100
#include <fstream>
std::ifstream fin("input.cnf");
std::ofstream fout("out35-17.txt");
double min, max, avg;
double foo[BITMAX];
int foon;
int dimension;
bool solutie[BITMAX];
bool populatie[1000][BITMAX];
double fitness[1000];
double rang[POPSIZE];
int popsize;
int totalLength;
bool newpop[POPSIZE][BITMAX];
void rand01();
void crossover();
void mutation();
void selection(int (*f)(bool* x, int n));
int select(double fs[], int size);
void geneticAlgorithm(int (*f)(bool* x, int n));

int cnf[BITMAX];
int cnfsize;
int fnc(bool* x, int n)
{
	int i;
	int result = 0;
	bool currentclause = false;
	for (i = 0; i < cnfsize; ++i)
	{
		if (cnf[i] > 0)
		{
			currentclause = currentclause || x[cnf[i] - 1];
		}
		else if (cnf[i] < 0)
		{
			currentclause = currentclause || !x[-cnf[i] - 1];
		}
		else
		{
			result += currentclause;
			currentclause = false;
		}
	}
	return result;
}
int main()
{
	max = MIN;
	int sample, i;
	dimension = 4;
	popsize = POPSIZE;
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << std::fixed << std::setprecision(5);
	char s[10];
	int clauze;
	fin >> s >> s >> s;
	fin >> dimension >> clauze;
	int x;
	while (!fin.eof())
	{
		fin >> x;
		cnf[cnfsize++] = x; 
	}
	for (sample = 0; sample < 30; ++sample)
	{
		foon = 0;
		geneticAlgorithm(fnc);
		std::cout << sample << ' ' << max << '\n';
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	for (int j = 0; j < totalLength; ++j)
		std::cout << solutie[j] << '\n';
	std::cout << std::fixed;
	fout << std::fixed;
	std::cout << "MIN: " << std::setprecision(5) << min << '\n';
	std::cout << "MAX: " << std::setprecision(5) << max << '\n';
	std::cout << "TIME: " << duration.count() << '\n';

	for (i = 0; i < foon; ++i)
	{
		avg += foo[i];
		fout << foo[i] << '\n';
	}
	avg /= (double)foon;
	double sd = 0;
	for (i = 0; i < foon; ++i)
	{
		sd += (foo[i] - avg) * (foo[i] - avg);
	}
	sd /= (double)(foon - 1);
	sd = sqrt(sd);
	std::cout << "Val No : " << foon << '\n';
	std::cout << "AVG: " << std::setprecision(5) << avg << '\n';
	std::cout << "SD: " << std::setprecision(5) << sd << '\n';
	return 0;
}

void rand01()
{
	int i, k;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> uniform_distance(0, 1);
	totalLength = dimension;
	for (k = 0; k < POPSIZE; ++k)
	{
		for (i = 0; i < totalLength; ++i)
		{
			populatie[k][i] = uniform_distance(mt);
		}
	}
}

void geneticAlgorithm(int (*f)(bool* x, int n))
{
	int t = 0;
	bool* xreal;
	int current_max;
	rand01();
	while (t < 500)
	{
		mutation();
		crossover();
		selection(f);
		int i;
		for (i = 0; i < popsize; ++i)
		{
			current_max = f(populatie[i], dimension);
			foo[foon++] = current_max;
			if (current_max > max)
			{
				max = current_max;
				for (int j = 0; j < totalLength; ++j)
					solutie[j] = populatie[i][j];
			}
			if (current_max < min)
				min = current_max;
		}
		t++;
	}
}
void crossover()
{
	bool* chromosome1 = new bool[BITMAX];
	bool* chromosome2 = new bool[BITMAX];
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> uniform_distance(1, totalLength - 1);
	std::uniform_real_distribution<double> rand01(0, 1);
	int pos = uniform_distance(mt);
	int i, k;

	for (i = 0; i < popsize; ++i)
	{
		rang[i] = rand01(mt);
	}

	bool sorting = true;
	double auxrang;
	bool* auxchromosome = new bool[BITMAX];
	while (sorting)
	{
		sorting = false;
		for (i = 0; i < popsize - 1; ++i)
		{
			if (rang[i] > rang[i + 1])
			{
				sorting = true;
				auxrang = rang[i + 1];
				for (k = 0; k < totalLength; ++k)
					auxchromosome[k] = populatie[i + 1][k];

				rang[i + 1] = rang[i];
				for (k = 0; k < totalLength; ++k)
					populatie[i + 1][k] = populatie[i][k];

				rang[i] = auxrang;
				for (k = 0; k < totalLength; ++k)
					populatie[i][k] = auxchromosome[k];
			}
		}
	}
	delete[] auxchromosome;

	bool mating = true;
	k = 0;
	while (mating)
	{
		if (rang[k] < 0.3 && rang[k + 1] < 0.3)
		{
			mating = true;
			k += 2;
			for (i = 0; i < pos; ++i)
			{
				chromosome1[i] = populatie[k][i];
				chromosome2[i] = populatie[k + 1][i];
			}
			for (i = pos; i < totalLength; ++i)
			{

				chromosome1[i] = populatie[k + 1][i];
				chromosome2[i] = populatie[k][i];
			}

			for (i = 0; i < totalLength; ++i)
			{
				populatie[popsize][i] = chromosome1[i];
				populatie[popsize + 1][i] = chromosome2[i];
			}
			popsize += 2;
		}
		else if (rang[k] < 0.3 && rang[k + 1] >= 0.3 && rand01(mt) < 0.5)
		{
			mating = false;
			k += 2;
			for (i = 0; i < pos; ++i)
			{
				chromosome1[i] = populatie[k][i];
				chromosome2[i] = populatie[k + 1][i];
			}
			for (i = pos; i < totalLength; ++i)
			{

				chromosome1[i] = populatie[k + 1][i];
				chromosome2[i] = populatie[k][i];
			}

			for (i = 0; i < totalLength; ++i)
			{
				populatie[popsize][i] = chromosome1[i];
				populatie[popsize + 1][i] = chromosome2[i];
			}
			popsize += 2;
		}
		else
		{
			mating = false;
		}
	}
	delete[] chromosome1;
	delete[] chromosome2;
}
void mutation()
{
	int k, i;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> uniform_distance(0, 1);
	for (k = 0; k < popsize; ++k)
	{
		for (i = 0; i < totalLength; ++i)
		{
			if (uniform_distance(mt) < 0.01)
				populatie[k][i] = !populatie[k][i];
		}
	}
}
void selection(int (*f)(bool* x, int n))
{
	int i, k, s;
	double epsilon = 0.000000001;
	for (i = 0; i < popsize; ++i)
	{
		fitness[i] = f(populatie[i], dimension) + epsilon;
	}

	double* fs = new double[popsize];


	fs[0] = fitness[0];

	for (i = 1; i < popsize; ++i)
	{
		fs[i] = fs[i - 1] + fitness[i];
	}

	for (i = 0; i < POPSIZE; ++i)
	{
		s = select(fs, popsize);
		for (k = 0; k < totalLength; ++k)
			newpop[i][k] = populatie[s][k];
	}
	popsize = POPSIZE;
	for (i = 0; i < POPSIZE; ++i)
	{
		for (k = 0; k < totalLength; ++k)
			populatie[i][k] = newpop[i][k];
	}
	delete[] fs;
}
int select(double fs[], int size)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> uniform_distance(0, 1);
	double pos = uniform_distance(mt) * fs[size - 1];

	for (int i = 0; i < size; ++i)
		if (pos <= fs[i])
			return i;
}