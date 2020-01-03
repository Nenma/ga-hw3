#include <iostream>
#include <math.h> 
#include <random>
#include <chrono> 
#include <iomanip>
#include <fstream>

#define MIN -9223372036854775
#define MAX 9223372036854775
#define BITMAX 100000

std::ifstream fin("input.cnf");
std::ofstream fout("out59-26.txt");
double min, max, avg;
double foo[BITMAX];
int foon;
double range[2][30];
int dimension;

int precision;
bool bitstring[BITMAX];
int length[30];
int totalLength;

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

bool solutie[BITMAX];
void annealing(int (*f)(bool* x, int n), double temp);
void rand01();

int main()
{
	min = MAX;
	max = MIN;
	int sample, i;

	char s[10];
	int clauze;
	fin >> s >> s >> s;
	fin >> dimension >> clauze;
	int x;
	cnfsize = 0;
	while (!fin.eof())
	{
		fin >> x;
		cnf[cnfsize++] = x;
	}
	totalLength = dimension;

	auto start = std::chrono::high_resolution_clock::now();
	for (sample = 0; sample < 30; ++sample)
	{
		annealing(fnc, 0.999);
		std::cout << sample << ' ' << max << '\n';
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
	std::cout << std::fixed;
	for (int i = 0; i < totalLength; ++i) {
		std::cout << solutie[i] << '\n';
	}
	std::cout << std::endl;
	std::cout << "MIN: " << std::setprecision(5) << min << '\n';
	std::cout << "MAX: " << std::setprecision(5) << max << '\n';
	std::cout << "TIME: " << duration.count() << '\n';
	fout << std::fixed;
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
	int i;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> uniform_distance(0, 1);
	for (i = 0; i < totalLength; ++i)
	{
		bitstring[i] = uniform_distance(mt);
	}
}

void annealing(int (*f)(bool* x, int n), double temp)
{
	rand01();
	double T, current_max;
	int i, pos;
	bool* b = new bool[BITMAX];
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> uniform_distance(0, totalLength - 1);
	std::uniform_real_distribution<double> random(0, 0.9999999);
	for (T = 100; T > 0.0001; T *= temp)
	{
		for (i = 0; i < totalLength; ++i)
		{
			b[i] = bitstring[i];
		}
		pos = uniform_distance(mt);
		bitstring[pos] = !bitstring[pos];

		if (f(b, dimension) > f(bitstring, dimension))
		{
			for (i = 0; i < totalLength; ++i)
			{
				bitstring[i] = b[i];
			}
		}
		else if (random(mt) > exp(-abs(f(b, dimension) - f(bitstring, dimension)) / T))
		{
			for (i = 0; i < totalLength; ++i)
			{
				bitstring[i] = b[i];
			}
		}
		current_max = f(bitstring, dimension);
		foo[foon++] = current_max;
		if (current_max > max) {
			max = current_max;
			for (int i = 0; i < totalLength; ++i) {
				solutie[i] = bitstring[i];
			}
		}
		if (current_max < min)
			min = current_max;
	}

	delete[] b;
}