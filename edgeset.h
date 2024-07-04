#ifndef _EDGESET
#define _EDGESET

#include<cassert>
#include<fstream>
#include<set>
using namespace std;

#define MAXSIZE 1000
#define MAXID 1000

typedef enum { root, divergence, reticulation, leaf } node_type;

ostream& operator<<(ostream &strm, const node_type &nt)
{
	switch (nt)
	{
		case 0: strm << "t"; break;
		case 1: strm << "d"; break;
		case 2: strm << "r"; break;
		case 3: strm << "l"; break;
	}

	return strm;
}

typedef enum { D, S, R, C, L, T } event_type;

ostream& operator<<(ostream &strm, const event_type &et)
{
	switch (et)
	{
		case 0: strm << "D"; break;
		case 1: strm << "S"; break;
		case 2: strm << "R"; break;
		case 3: strm << "C"; break;
		case 4: strm << "L"; break;
		case 5: strm << "T"; break;
	}

	return strm;
}


typedef struct node node;

struct node
{
	int id;
	string name;
	double date, breakpoint;
	node_type type;
	node *parent[2];
	node *child[2];
} ;

typedef struct 
{
	node *f[MAXID];
	event_type e[MAXID];
} reconciliation;


/******/

class edge
{
	public:
		node* parent;
		node* child;
		edge();
		edge(node*, node*);
} ;

edge::edge()
{
	parent = child = NULL;
}

edge::edge(node *p, node *c)
{
	parent = p;
	child = c;
}

bool operator<(const edge &e1, const edge &e2)
{
	if (e1.parent == NULL)
		return false;
	
	if (e2.parent == NULL)
		return true;

	if (e1.child->id < e2.child->id)
		return true;
	
	if (e1.child->id > e2.child->id)
		return false;

	if (e1.parent->id < e2.parent->id)
		return true;

	return false;
}

bool operator==(const edge &e1, const edge &e2)
{
	//Assumes no double edges!
	return (e1.child == e2.child && e1.parent == e2.parent);
}

ostream& operator<<(ostream &strm, const edge &e)
{
	strm << "(";
	if (e.parent == NULL)
		strm << "x";
	else
		strm << e.parent->id;
	strm << ", " << e.child->id << ")";
	return strm;
}

typedef struct map map;


/******/

node* read_gene_network(char* file, node** genes, int &ngenes)
//Read gene network from file
//Returns root of gene network
//Still assumes genes are in bottom-up order :( also that gene id corresponds to position
{
	ifstream in(file);

	int nt, id;
	bool has_dates;

	in >> ngenes >> has_dates;

	for (int i = 0; i < ngenes; i++)
		genes[i] = new node;

	for (int i = 0; i < ngenes; i++)
	{
		in >> genes[i]->id >> nt;
		genes[i]->type = (node_type)nt;

		if (genes[i]->type != root)
		{
			in >> id;
			genes[i]->parent[0] = genes[id];

			if (genes[i]->type == reticulation)
			{
				in >> id;
				genes[i]->parent[1] = genes[id];
				in >> genes[i]->breakpoint;
			}
			else
				genes[i]->parent[1] = NULL;
		}

		if (genes[i]->type == leaf)
			in >> genes[i]->name;

		if (has_dates)
			in >> genes[i]->date;
		else
			genes[i]->date = 0;
	}

	//Set children
	for (int i = 0; i < ngenes; i++)
	{
		if (genes[i]->type == root)
			genes[i]->parent[0] = genes[i]->parent[1] = NULL;
		else
		{
			if (genes[genes[i]->parent[0]->id]->child[0] == NULL)
				genes[genes[i]->parent[0]->id]->child[0] = genes[i];
			else
				genes[genes[i]->parent[0]->id]->child[1] = genes[i];

			if (genes[i]->type == reticulation)
			{
				if (genes[genes[i]->parent[1]->id]->child[0] == NULL)
					genes[genes[i]->parent[1]->id]->child[0] = genes[i];
				else
					genes[genes[i]->parent[1]->id]->child[1] = genes[i];
			}
		}
	}

	in.close();

	return genes[ngenes-1];
}


node* read_species_tree(char* file, node** species, int &nspecies)
{
	node* s = read_gene_network(file, species, nspecies);

	if (s->date == 0)
	{
		//Set dates automatically
		for (int i = 0; i < nspecies; i++)
		{
			if (species[i]->type == leaf)
				species[i]->date = 0;
			else if (species[i]->type == divergence)
				species[i]->date = max(species[i]->child[0]->date, species[i]->child[1]->date) + 1;
			else
				species[i]->date = species[i]->child[0]->date + 1;
		}
	}

	return s;
}

void assign_dates_helper(node* n, set<node*> &seen)
{
	if (seen.find(n) != seen.end())
		return;

	switch (n->type)
	{
		case leaf: 
			n->date = 0;
			break;
		case root: case reticulation:
			assign_dates_helper(n->child[0], seen);
			n->date = n->child[0]->date + 1;
			break;
		case divergence:
			assign_dates_helper(n->child[0], seen);
			assign_dates_helper(n->child[1], seen);
			n->date = max(n->child[0]->date, n->child[1]->date) + 1;
			break;
	}

	seen.insert(n);
}

void assign_dates(node* n)
	//assigns arbitrary but self-consistent (parent > child) dates
{
	set<node*> seen;
	seen.clear();
	assign_dates_helper(n, seen);
}

void get_leaves_helper(node* n, set<node*> &leaves, set<node*> &seen)
{
	if (seen.find(n) != seen.end())
		return;

	switch (n->type)
	{
		case leaf: 
			if (leaves.find(n) == leaves.end())
				leaves.insert(n); break;
		case root: case reticulation:
			get_leaves_helper(n->child[0], leaves, seen);
			break;
		case divergence:
			get_leaves_helper(n->child[0], leaves, seen);
			get_leaves_helper(n->child[1], leaves, seen);
			break;
	}

	seen.insert(n);
}

void get_leaves(node* n, set<node*> &leaves)
{
	set<node*> seen;
	seen.clear();
	get_leaves_helper(n, leaves, seen);
}

double find_paralogy(node* l, node* m, const reconciliation &r, double lhs)
	//assume network is dated with assign_dates
{
	double rhs = 1;
	
	node *parent_l = l, *parent_m = m;
	node *increase;

	while (parent_l != parent_m)
	{
		if (parent_l->date < parent_m->date)
		{
			if (parent_l->type == reticulation)
			{
				if (parent_l->breakpoint <= lhs)
					parent_l = parent_l->parent[1];
				else
				{
					rhs = min(rhs, parent_l->breakpoint);
					parent_l = parent_l->parent[0];
				}
			}
			else
				parent_l = parent_l->parent[0];
		}
		else
		{
			if (parent_m->type == reticulation)
			{
				if (parent_m->breakpoint <= lhs)
					parent_m = parent_m->parent[1];
				else
				{
					rhs = min(rhs, parent_m->breakpoint);
					parent_m = parent_m->parent[0];
				}
			}
			else
				parent_m = parent_m->parent[0];
		}
	}

	cout << l->id << '\t' << m->id << '\t' << lhs << '\t' << rhs << '\t' << (r.e[parent_l->id] == S ? "S" : "D") << endl;

	return rhs;
}

void print_paralogy_info(node* n, const reconciliation &r)
{
	set<node*> leaves;

	assign_dates(n);
	get_leaves(n, leaves);

	cout << "Gene paralogy information" << endl;

	for (set<node*>::iterator l = leaves.begin(); l != leaves.end(); l++)
	{
		set<node*>::iterator m = l;
		m++;

		for (; m != leaves.end(); m++)
		{
			double rhs = find_paralogy(*l, *m, r, 0);

			while (rhs < 1)
				rhs = find_paralogy(*l, *m, r, rhs);
		}
	}

	cout << endl;
}

/******/

class edgeset
{
	public:
		edge edges[MAXSIZE];
		int nedges;
		edgeset();
		void add(edge);
		void expand(int);			//Expand edge to reticulation
		void merge(int, int);		//Merge duplication
		void operator=(const edgeset&);

} ;

edgeset::edgeset()
{
	nedges = 0;
}

void edgeset::add(edge e)
{
	edges[nedges].child = e.child;
	edges[nedges].parent = e.parent;
	nedges++;
}

void edgeset::expand(int i)
{
	assert(i < nedges && edges[i].parent->type == reticulation);

	edges[nedges].child = edges[i].parent;
	edges[nedges].parent = edges[i].parent->parent[1];
	edges[i].child = edges[i].parent;
	edges[i].parent = edges[i].parent->parent[0];
	nedges++;
}

void edgeset::merge(int i, int j)
{
	assert(i < nedges && j < nedges && edges[i].parent == edges[j].parent && (edges[i].parent->type == divergence || edges[i].parent->type == root));

	edges[i].child = edges[i].parent;
	edges[i].parent = edges[i].parent->parent[0];
	edges[j].child = edges[nedges-1].child;
	edges[j].parent = edges[nedges-1].parent;
	nedges--;
}

void edgeset::operator=(const edgeset &src)
{
	nedges = src.nedges;

	for (int i = 0; i < nedges; i++)
	{
		edges[i].parent = src.edges[i].parent;
		edges[i].child = src.edges[i].child;
	}
}

edgeset operator+(const edgeset &s1, const edgeset &s2)
{
	edgeset result;

	result = s1;

	bool merged;

	for (int i = 0; i < s2.nedges; i++)
	{
		merged = false;

		for (int j = 0; j < s1.nedges; j++)
			if (result.edges[j].parent == s2.edges[i].parent)
			{
				//Merge
				result.edges[j].child = result.edges[j].parent;
				result.edges[j].parent = result.edges[j].parent->parent[0];

				merged = true;

				break;
			}

		if (!merged)
		{
			//Add
			result.edges[result.nedges].child = s2.edges[i].child;
			result.edges[result.nedges].parent = s2.edges[i].parent;
			result.nedges++;
		}
	}

	return result;
}


#endif
