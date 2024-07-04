#include<iostream>
#include<random>
#include<cstdlib>
#include<cfloat>
#include<set>
#include<cstring>
using namespace std;

#ifndef MAXID
#define MAXID 1000
#endif

double Drate = 0.5;
double Lrate = 0.5;
double Rrate = 0.5;

#include"edgeset.h"

void find_next_event(node** current, int num_cur, const reconciliation &recon, double time, int &next_select, int &rec_select, event_type& next_event, double& next_time, double& rec_breakpoint);
bool apply_loss(node* n);
void relabel_ids(node* n, set<node*> &seen, reconciliation &r, const reconciliation &recon, int &num_id);
void print_network(node* n, set<node*> &seen);
void print_reconciliation(node* n, set<node*> &seen, const reconciliation &r);

int main(int argc, char** argv)
{
	node *s;
	node *species[MAXID];
	int nspecies;
	s = read_species_tree(argv[1], species, nspecies);

	//set rates from command line options
	for (int i = 2; i < argc; i += 2)
	{
		if (strcmp(argv[i],"--Drate") == 0)
			Drate = atof(argv[i+1]);
		else if (strcmp(argv[i],"--Lrate") == 0)
			Lrate = atof(argv[i+1]);
		else if (strcmp(argv[i],"--Rrate") == 0)
			Rrate = atof(argv[i+1]);
	}

	node *rt;

	//Generate gene network
	reconciliation recon;
	node* current[MAXID];
	int num_id = 0, num_cur = 0;
	double time = s->date;

	rt = new node;
	recon.f[num_id] = s->child[0];
	recon.e[num_id] = T;
	rt->id = num_id++;
	rt->date = time;
	rt->type = root;

	rt->child[0] = new node;
	recon.f[num_id] = s->child[0];
	current[num_cur++] = rt->child[0];
	rt->child[0]->parent[0] = rt;
	rt->child[0]->id = num_id++;

	event_type next_event;
	double next_time, rec_breakpoint;
	int next_select, rec_select;
	node* next_gene, *rec_gene;

	node* next_species;
	int old_cur;

	//Find next event
	find_next_event(current, num_cur, recon, time, next_select, rec_select, next_event, next_time, rec_breakpoint);

	while (next_event != C)
	{
		time = next_time;
		next_gene = current[next_select];

		switch(next_event)
		{
			case S:
				next_species = recon.f[next_gene->id];
				old_cur = num_cur;

				for (int i = 0; i < old_cur; i++)
					if (recon.f[current[i]->id] == next_species)
					{
						current[i]->type = divergence;
						current[i]->date = next_time;
						recon.e[current[i]->id] = S;

						for (int j = 0; j < 2; j++)
						{
							current[i]->child[j] = new node;

							recon.f[num_id] = recon.f[current[i]->id]->child[j];
							current[i]->child[j]->id = num_id++;
							current[i]->child[j]->parent[0] = current[i];
						}

						current[num_cur++] = current[i]->child[1];
						current[i] = current[i]->child[0];
					}
				break;
			case D:
				next_gene->type = divergence;
				next_gene->date = next_time;
				recon.e[next_gene->id] = D;

				for (int i = 0; i < 2; i++)
				{
					next_gene->child[i] = new node;

					recon.f[num_id] = recon.f[next_gene->id];

					next_gene->child[i]->id = num_id++;
					next_gene->child[i]->parent[0] = next_gene;
				}

				current[next_select] = next_gene->child[0];
				current[num_cur++] = next_gene->child[1];
				break;

			case R:
				//Need to accept that gene ids are not dense
				//But <current> should be
				rec_gene = current[rec_select];
				//cout << next_gene->id << "\t" << rec_gene->id << endl;

				if (next_gene->parent[0] == rec_gene->parent[0])
				{
					//Stitching an eyelet - remove? was removed but I put it back in (14/3/24)
					recon.f[next_gene->parent[0]->id] = recon.f[next_gene->id];
					current[next_select] = next_gene->parent[0];
					current[rec_select] = current[--num_cur];
					delete next_gene;
					delete rec_gene;

					break;
				}

				next_gene->type = reticulation;
				next_gene->date = next_time;
				next_gene->breakpoint = rec_breakpoint;
				recon.e[next_gene->id] = R;

				next_gene->parent[1] = rec_gene->parent[0];
				rec_gene->parent[0]->child[rec_gene->parent[0]->child[0] == rec_gene ? 0 : 1] = next_gene;
				next_gene->child[0] = new node;

				recon.f[num_id] = recon.f[next_gene->id];

				next_gene->child[0]->id = num_id++;
				next_gene->child[0]->parent[0] = next_gene;

				current[next_select] = next_gene->child[0];
				current[rec_select] = current[--num_cur];

				delete rec_gene;
				break;
			case L:
				bool eyelet = apply_loss(next_gene);
				current[next_select] = current[--num_cur];

				if (num_cur == 0 || eyelet)
				{
					//Last gene lost, start again
					cout << "Last gene lost or eyelet, restarting" << endl;
					num_id = 0;
					num_cur = 0;
					time = s->date;

					rt = new node;
					recon.f[num_id] = s;
					rt->id = num_id++;
					rt->date = time;
					rt->type = root;

					rt->child[0] = new node;
					recon.f[num_id] = s->child[0];
					current[num_cur++] = rt->child[0];
					rt->child[0]->parent[0] = rt;
					rt->child[0]->id = num_id++;
				}

				break;
		}

		//Find next event
		find_next_event(current, num_cur, recon, time, next_select, rec_select, next_event, next_time, rec_breakpoint);
	}

	//Set leaves
	for (int i = 0; i < num_cur; i++)
	{
		current[i]->type = leaf;
		current[i]->child[0] = current[i]->child[1] = NULL;
		current[i]->name = recon.f[current[i]->id]->name;
		recon.e[current[i]->id] = C;
	}

	//Need to relabel ids and then output
	reconciliation new_recon;

	set<node*> seen;
	seen.clear();

	int new_id = 0;
	relabel_ids(rt, seen, new_recon, recon, new_id);
	cout << endl << "Gene network" << endl;
	cout << new_id << " 0" << endl;
	seen.clear();
	print_network(rt, seen);
	cout << endl;

	//Print reconciliation
	seen.clear();
	cout << "Reconciliation" << endl;
	print_reconciliation(rt, seen, new_recon);
	cout << endl;

	//Print paralogy information
	print_paralogy_info(rt, new_recon);

	return 0;
}

void find_next_event(node** current, int num_cur, const reconciliation &recon, double time, int &next_select, int &rec_select, event_type& next_event, double& next_time, double& rec_breakpoint)
{
//	cout << "current: " << endl;
//	for (int i = 0; i < num_cur; i++)
//		cout << current[i]->id << " " << recon.f[current[i]->id]->id << " | ";
//	cout << endl;

	//Finds next event type, time and gene
	random_device rd;
	mt19937 gen(rd());
	//mt19937 gen(seed++);
	
	//Calculate next speciation time and lineage
	double next_S = 0;
	int next_S_select = 0;

	for (int i = 0; i < num_cur; i++)
	{
		if (recon.f[current[i]->id]->date > next_S)
		{
			next_S = recon.f[current[i]->id]->date;
			next_S_select = i;
		}
	}

	//Calculate R rate
	//Currently rate is Rrate for each gene AFTER the first in a species
	int total_R = 0;
	int species_R[MAXID];

	for (int i = 0; i < MAXID; i++)
		species_R[i] = -1;

	for (int i = 0; i < num_cur; i++)
		species_R[recon.f[current[i]->id]->id]++;

	for (int i = 0; i < MAXID; i++)
		total_R += max(0, species_R[i]);

	//Calculate D/L/R rate
	double DLRrate = num_cur * (Drate + Lrate) + total_R * Rrate;

	//Calculate next time
	exponential_distribution<> exp(DLRrate);

	next_time = time - exp(gen);

	if (next_time < next_S)
	{
		//Speciation or leaf
		if (recon.f[current[next_S_select]->id]->type == leaf)
		{
			//Leaf
			next_event = C;
			next_time = 0;
		}
		else
		{
			//Speciation
			next_event = S;
			next_time = next_S;
			next_select = next_S_select;
		}

//		cout << "next event is " << next_event << " at time " << next_time << endl;
		return;
	}

	//Otherwise, D/L/R
	//Select event type
	uniform_real_distribution<> unif_real(0, DLRrate);
	double select = unif_real(gen);

	if (select < num_cur*Drate)
		next_event = D;
	else if (select < num_cur*(Drate + Lrate))
		next_event = L;
	else
		next_event = R;

//	cout << "next event is " << next_event << " at time " << next_time << endl;

	//Select event gene
	if (next_event != R)
	{
		uniform_int_distribution<> unif(0, num_cur-1);
		next_select = unif(gen);
	}
	else	
	{
		uniform_int_distribution<> unif(1, total_R);
		select = unif(gen);
		int sp_id = 0;

		for (sp_id = 0; sp_id < MAXID && select > 0; sp_id++)
			select -= max(0,species_R[sp_id]);
		sp_id--;

		uniform_int_distribution<> unif2(1, species_R[sp_id]+1);
		next_select = unif2(gen);

		int gene_id;
		for (gene_id = 0; gene_id < num_cur && next_select > 0; gene_id++)
			if (recon.f[current[gene_id]->id]->id == sp_id)
				next_select--;

		next_select = gene_id-1;

		uniform_int_distribution<> unif3(1, species_R[sp_id]);
		rec_select = unif3(gen);

		for (gene_id = 0; gene_id < num_cur && rec_select > 0; gene_id++)
			if (recon.f[current[gene_id]->id]->id == sp_id && gene_id != next_select)
				rec_select--;

		rec_select = gene_id-1;

		uniform_real_distribution<> unif4(0.0, 1.0);
		rec_breakpoint = unif4(gen);
	}
}

bool apply_loss(node* n)
{
	//By convention, a loss can only be applied to a "current" node, which has no type and only 1 parent (in index 0)
	bool eyelet = false;

	node* sibling;
	if (n->parent[0]->type == divergence)
		sibling = n->parent[0]->child[n->parent[0]->child[0] == n ? 1 : 0];

	switch(n->parent[0]->type)
	{
		case root:
			delete n->parent[0];
			break;
		case divergence:
			sibling->parent[sibling->parent[0] == n->parent[0] ? 0 : 1] = n->parent[0]->parent[0];
			n->parent[0]->parent[0]->child[n->parent[0]->parent[0]->child[0] == n->parent[0] ? 0 : 1] = sibling;

			//Some code for if this results in an eyelet (reticulation node with two parents the same)
			if (sibling->type == reticulation && sibling->parent[0] == sibling->parent[1])
				eyelet = true;
				//sibling->child[0]->parent[sibling->child[0]->parent[0] == sibling ? 0 : 1] = sibling->parent[0]->parent[0];

			delete n->parent[0];
			break;
		case reticulation:
			node* temp = new node;
			temp->parent[0] = n->parent[0]->parent[1];
			n->parent[0]->parent[1]->child[n->parent[0]->parent[1]->child[0] == n->parent[0] ? 0 : 1] = temp;
			apply_loss(n->parent[0]);
			apply_loss(temp);
			break;
	}

	delete n;

	return eyelet;
}

//Not a tree so need to keep list of seen nodes!
//Are we going to have a problem because the id of the nodes change and therefore the comparison operator may change? Not sure, but keep in mind
void relabel_ids(node* n, set<node*> &seen, reconciliation &r, const reconciliation &recon, int &num_id)
{
	if (seen.find(n) != seen.end())
		return;

	//Post hoc eyelet removal - apprently does not completely work?
	if (n->type == divergence && n->child[0] == n->child[1])
	{
		n->parent[0]->child[n->parent[0]->child[0] == n ? 0 : 1] = n->child[0]->child[0];
		n->child[0]->child[0]->parent[n->child[0]->child[0]->parent[0] == n->child[0] ? 0 : 1] = n->parent[0];
		relabel_ids(n->child[0]->child[0], seen, r, recon, num_id);
		return;
	}

	//Relabel ids in bottom-up order (and dense)
	//Also relabel species_map
	switch (n->type)
	{
		case root: case reticulation:
			relabel_ids(n->child[0], seen, r, recon, num_id);
			break;
		case divergence:
			relabel_ids(n->child[0], seen, r, recon, num_id);
			relabel_ids(n->child[1], seen, r, recon, num_id);
			break;
	}

	r.f[num_id] = recon.f[n->id];
	r.e[num_id] = recon.e[n->id];
	n->id = num_id++;
	seen.insert(n);
}

void print_network(node* n, set<node*> &seen)
{
	if (seen.find(n) != seen.end())
		return;

	//Print out network in bottom-up order
	switch (n->type)
	{
		case root: case reticulation:
			print_network(n->child[0], seen);
			break;
		case divergence:
			print_network(n->child[0], seen);
			print_network(n->child[1], seen);
			break;
	}

	cout << n->id << " " << (int)n->type;
	switch(n->type)
	{
		case divergence:
			cout << " " << n->parent[0]->id;
			break;
		case reticulation:
			cout << " " << n->parent[0]->id << " " << n->parent[1]->id << " " << n->breakpoint;
			break;
		case leaf: 
			cout << " " << n->parent[0]->id << " " << n->name;
			break;
	}
	cout << endl;
	
	seen.insert(n);
}

void print_reconciliation(node* n, set<node*> &seen, const reconciliation &r)
{
	if (seen.find(n) != seen.end())
		return;

	switch (n->type)
	{
		case root: case reticulation:
			print_reconciliation(n->child[0], seen, r);
			break;
		case divergence:
			print_reconciliation(n->child[0], seen, r);
			print_reconciliation(n->child[1], seen, r);
			break;
	}

	cout << n->id << " " << r.f[n->id]->id << " " << r.e[n->id] << endl;

	seen.insert(n);
}

/* Need to output true reconciliation								./
 * Obviously as well as gene network (in separate file)
 * Did I stuff up time allocation in the species tree? Yes, but it's fixed	./
 */

//Problem: can output network with double edges (assumed not to happen in perecon)
//10/4/24: this is now fixed (I hope)
