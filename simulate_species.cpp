#include<iostream>
#include<random>
#include<cstdlib>
using namespace std;

#ifndef MAXID
#define MAXID 1000
#endif

#define Srate 1

#include"edgeset.h"


int main(int argc, char** argv)
{
	int number_of_species = (int) strtol(argv[1], (char **)NULL, 10);

	random_device rd;
	mt19937 gen(rd());

	double time = 0;

	node *rt, *child[2];

	//Generate species tree
	node* current[MAXID];
	node* all[MAXID];
	int num_id = 0, num_cur = 0;

	rt = new node;
	all[num_id] = rt;
	rt->id = num_id++;
	rt->date = time;
	rt->type = root;

	rt->child[0] = new node;
	all[num_id] = rt->child[0];
	current[num_cur++] = rt->child[0];
	rt->child[0]->parent[0] = rt;
	rt->child[0]->id = num_id++;

	//Branching process
	while (num_cur < number_of_species)
	{
		exponential_distribution<> exp(num_cur*Srate);
		uniform_int_distribution<> unif(0,num_cur-1);

		time += exp(gen);

		int select = unif(gen);

		for (int i = 0; i < 2; i++)
		{
			child[i] = new node;
			child[i]->id = num_id++;

			current[select]->child[i] = child[i];
			child[i]->parent[0] = current[select];

			all[num_id-1] = child[i];
		}

		current[select]->date = time;
		current[select]->type = divergence;
		current[select] = child[0];
		current[num_cur++] = child[1];
	}

	exponential_distribution<> exp_last(num_cur*Srate);
	time += exp_last(gen);

	//Set leaves
	for (int i = 0; i < num_cur; i++)
	{
		current[i]->type = leaf;
		current[i]->child[0] = current[i]->child[1] = NULL;
		current[i]->date = time;
	}

	//Change ids to bottom-up order - id goes to num_id - id - 1

	//Output, in bottom-up order
	cout << num_id << " 1" << endl;
	char name = 'a';
	for (int i = num_id-1; i >= 0; i--)
	{
		cout << num_id - all[i]->id - 1 << " " << (int)all[i]->type;
		if (all[i]->type != root)
			cout << " " << num_id - all[i]->parent[0]->id - 1;
		if (all[i]->type == leaf)
		{
			all[i]->name = name++;
			cout << " " << all[i]->name;
		}
		cout << " " << (all[num_id-1]->date - all[i]->date) << endl;
	}

	return 0;
}

/* Still need to get times in the file somehow		./
 * Include simulation of genes - different file
 */
