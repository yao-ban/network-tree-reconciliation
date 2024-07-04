#include<iostream>
#include<fstream>
#include<cstdlib>
#include<climits>
#include<set>
using namespace std;

#include"edgeset.h"

#define MAXID 1000
#define Lcost 1
#define Dcost 1

node* create_species_tree(node** species, node* parent, int level, int &id);
int dist(node* s1, node* s2);
void calculate_lca(node** genes, int ngenes, node** species_map, reconciliation &r);
void calculate_lca_hca(node** genes, int ngenes, node** species_map, reconciliation &r);
void handle_lca_hca(int i, node** genes, reconciliation &lca, reconciliation &r);
void decompose_bc(node** genes, int ngenes, node** rbc, int* bc_num, bool* is_treechild, int &nbc);
void decompose_bc_helper(node** genes, int ngenes, bool* is_rbc, int* depth, int* low, node** parent, node* g, int d);
void calculate_mpnsr(node** genes, int ngenes, node** species, int nspecies, node** species_map, reconciliation &r);
void calculate_mpr(node** genes, int ngenes, node** species, int nspecies, node** species_map, reconciliation &r);
void resolve_backtracking(map &top, node *species, int bcc, int* bc_num, reconciliation &r);
void set_from(node *gene, map **back, node *species, int *bc_num, reconciliation &r);
void create_edge_list(node** genes, int ngenes, edge** edges, int &nedges);

int main(int argc, char** argv)
{
	node *g, *s;

	node *genes[MAXID], *species[MAXID], *species_map[MAXID];

	int ngenes = 0;
	g = read_gene_network(argv[1], genes, ngenes);

	int nspecies = 0;
	s = read_species_tree(argv[2], species, nspecies);

	cout << "species tree" << endl;
	for (int i = 0; i < nspecies; i++)
		cout << species[i]->id << '\t' << species[i]->type << '\t' << ( species[i]->type == root ? -1 : species[i]->parent[0]->id ) << '\t' << ( species[i]->type == leaf ? -1 : species[i]->child[0]->id ) << '\t' << ( species[i]->type == divergence ? species[i]->child[1]->id : -1 ) << endl;
	cout << endl;

	//Automatic species mapping
	for (int i = 0; i < ngenes; i++)
		if (genes[i]->type == leaf)
		{
			for (int j = 0; j < nspecies; j++)
				if (genes[i]->name == species[j]->name)
				{
					species_map[i] = species[j];
					break;
				}
		}

	cout << "gene network" << endl;
	cout << "id\ttype\tp(s)\tc(s)\tspecies" << endl;
	for (int i = 0; i < ngenes; i++)
	{
		if (genes[i]->type != reticulation)
			cout << genes[i]->id << '\t' << genes[i]->type << '\t' << ( genes[i]->type == root ? -1 : genes[i]->parent[0]->id ) << '\t' << ( genes[i]->type == leaf ? -1 : genes[i]->child[0]->id ) << '\t' << ( genes[i]->type == divergence ? genes[i]->child[1]->id : -1 ) << '\t' << ( genes[i]->type == leaf ? species_map[i]->id : -1 ) << endl;
		else
			cout << genes[i]->id << '\t' << genes[i]->type << '\t' << genes[i]->parent[0]->id << '\t' << genes[i]->parent[1]->id << '\t' << genes[i]->child[0]->id << '\t' << -1 << endl;
	}
	cout << endl;

	/******/

	reconciliation lca, lca_hca;

	calculate_lca(genes, ngenes, species_map, lca);

	cout << "LCA reconciliation" << endl;
	for (int i = 0; i < ngenes; i++)
		cout << genes[i]->id << '\t' << lca.f[i]->id << '\t' << lca.e[i] << endl;
	cout << endl;

	calculate_lca_hca(genes, ngenes, species_map, lca_hca);

	cout << "LCA-HCA reconciliation" << endl;
	for (int i = 0; i < ngenes; i++)
		cout << genes[i]->id << '\t' << lca_hca.f[i]->id << '\t' << lca_hca.e[i] << endl;
	cout << endl;

	/******/

	reconciliation mpr;

	calculate_mpr(genes, ngenes, species, nspecies, species_map, mpr);

	cout << "MPR" << endl;
	for (int i = 0; i < ngenes; i++)
		cout << genes[i]->id << '\t' << mpr.f[i]->id << '\t' << mpr.e[i] << endl;
	cout << endl;
	
	print_paralogy_info(g, mpr);

	return 0;

}

node* create_species_tree(node** species, node* parent, int level, int &id)
//Create species tree. Currently hard-coded to balanced tree with <level> levels
//Returns root of species tree
//No longer actually use this function but I'll leave it for legacy purposes - may need some changes though
{
	node *n = new node;

	n->parent[0] = parent;

	if (level == 0)
		n->type = leaf;
	else
	{
		for (int i = 0; i < 2; i++)
			n->child[i] = create_species_tree(species, n, level-1, id);

		n->type = divergence;
	}

	n->id = id;
	n->date = level;
	species[id] = n;
	id++;

	return n;
}

int dist(node* s1, node* s2)
	//Will fail (hang) if s1 and s2 are not comparable
{
	int dist = 0;
	node *s, *end;

	if (s1->id < s2->id)
	{
		s = s1;
		end = s2;
	}
	else
	{
		s = s2;
		end = s1;
	}

	while (s != end)
	{
		s = s->parent[0];
		dist++;
	}

	return dist;
}

void calculate_lca(node** genes, int ngenes, node** species_map, reconciliation &r)
//Calculate LCA reconciliation
//Reconciliation is an array of mappings indexed by gene id(!) - not necessarily the same order as the gene array, although currently it is
//Requires bottom-up order in genes array
//Requires date in species tree, but only requires date consistency between parent/child
{
	node* child_lca[2];

	for (int i = 0; i < ngenes; i++)
	{
		switch (genes[i]->type)
		{
			case leaf: 	
				r.f[i] = species_map[i];
				r.e[i] = C;
				break;

			case reticulation:
				r.f[i] = r.f[genes[i]->child[0]->id];
				r.e[i] = R;
				break;

			case divergence:
				for (int j = 0; j < 2; j++)
					child_lca[j] = r.f[genes[i]->child[j]->id];

				while (child_lca[0] != child_lca[1])
				{
					if (child_lca[0]->date < child_lca[1]->date)
						child_lca[0] = child_lca[0]->parent[0];
					else
						child_lca[1] = child_lca[1]->parent[0];
				}

				r.f[i] = child_lca[0];
				if (r.f[i] == r.f[genes[i]->child[0]->id] || r.f[i] == r.f[genes[i]->child[1]->id])
					r.e[i] = D;
				else
					r.e[i] = S;
				break;

			case root:
				r.f[i] = r.f[genes[i]->child[0]->id];
				r.e[i] = T;
				break;
		}
	}
}

void calculate_lca_hca(node** genes, int ngenes, node** species_map, reconciliation &r)
//Calculate LCA-HCA reconciliation
{
	reconciliation lca;

	calculate_lca(genes, ngenes, species_map, lca);

	for (int i = 0; i < ngenes; i++)
		r.f[i] = NULL;

	for (int i = 0; i < ngenes; i++)
		if (genes[i]->type != reticulation)
		{
			r.f[i] = lca.f[i];
			r.e[i] = lca.e[i];
		}

	for (int i = 0; i < ngenes; i++)
		if (r.f[i] == NULL)
			handle_lca_hca(i, genes, lca, r);
}

void handle_lca_hca(int i, node** genes, reconciliation &lca, reconciliation &r)
//Subsidiary function for reticulation nodes for LCA-HCA reconciliation
{
	//Calculate LCA-HCA of parents
	if (r.f[genes[i]->parent[0]->id] == NULL)
		handle_lca_hca(genes[i]->parent[0]->id, genes, lca, r);
	if (r.f[genes[i]->parent[1]->id] == NULL)
		handle_lca_hca(genes[i]->parent[1]->id, genes, lca, r);

	node* parent_lca[2];
	node *node_lca, *prev_node;
	for (int j = 0; j < 2; j++)
	{
		node_lca = lca.f[i];
		parent_lca[j] = r.f[genes[i]->parent[j]->id];

		if (r.e[genes[i]->parent[j]->id] == S)
		{
			//move down - have to decide which child to move to
			while (node_lca != parent_lca[j])
			{
				prev_node = node_lca;
				node_lca = node_lca->parent[0];
			}

			parent_lca[j] = prev_node;
		}
	}

	r.f[i] = ( parent_lca[0]->date < parent_lca[1]->date ? parent_lca[0] : parent_lca[1] );
	r.e[i] = R;
}

void decompose_bc(node** genes, int ngenes, node** rbc, int* bc_num, bool* is_treechild, int &nbc)
//Calculate the root of the biconnected component for each node
//Store as pointer in <rbc>
//<bc_num> stores label for each non-leaf bcc
//<is_treechild> determines if each bcc is tree-child or not, indexed by bcc num
//Still assumes that gene id is equal to position in gene array :(
{
	int depth[MAXID], low[MAXID];
	bool is_rbc[MAXID];
	node* parent[MAXID];

	for (int i = 0; i < ngenes; i++)
	{
		is_rbc[i] = (genes[i]->type == root ? true : false);
		low[i] = MAXID+1;
		depth[i] = -1;
	}

	//Find articulation points
	decompose_bc_helper(genes, ngenes, is_rbc, depth, low, parent, genes[0], 0);

	//Find roots of biconnected components from top down
	//Also label biconnected components which are not leaves, ordered bottom-up
	nbc = 0;
	for (int i = 0; i < ngenes; i++)
		bc_num[i] = 1;

	for (int i = ngenes-1; i >= 0; i--)
	{
		if (is_rbc[i])
		{
			rbc[i] = genes[i];
			if (genes[i]->type != leaf)
				bc_num[i] = -(nbc++);
		}
		else
		{
			rbc[i] = rbc[genes[i]->parent[0]->id];
			bc_num[i] = bc_num[genes[i]->parent[0]->id];
		}
	}

	for (int i = 0; i < ngenes; i++)
	{
		if (bc_num[i] == 1)
			bc_num[i] = -1;
		else
			bc_num[i] += nbc-1;
	}

	//Determine if each bcc is tree-child
	for (int b = 0; b < nbc; b++)
		is_treechild[b] = true;

	for (int i = 0; i < ngenes; i++)
		if (genes[i]->type != leaf)
		{
			if (genes[i]->type == reticulation && genes[i]->child[0]->type == reticulation)
				is_treechild[bc_num[i]] = false;
			else if (genes[i]->child[0]->type == reticulation && genes[i]->child[1]->type == reticulation)
				is_treechild[bc_num[i]] = false;
		}

	cout << "BCC decomposition" << endl;
	cout << nbc << " non-leaf biconnected components" << endl;
	for (int i = 0; i < ngenes; i++)
		cout << genes[i]->id << '\t' << depth[i] << '\t' << low[i] << '\t' << (is_rbc[i] ? "t" : "f") << '\t' << rbc[i]->id << '\t' << bc_num[i] << endl;
	cout << endl;
	cout << "Are BCCs tree-child?" << endl;
	for (int i = 0; i < nbc; i++)
		cout << i << '\t' << (is_treechild[i] ? "t" : "f") << endl;
	cout << endl;

}

void decompose_bc_helper(node** genes, int ngenes, bool* is_rbc, int* depth, int* low, node** parent, node* g, int d)
{
	depth[g->id] = d;
	low[g->id] = d;

	node* children[4];
	int nc = 0;
	bool seen_parent = false;

	//<children> is the children of the node <g> in the DFS search tree (NOT the original tree/network)
	if (g->type != root)
		children[nc++] = g->parent[0];
	if (g->type == reticulation)
		children[nc++] = g->parent[1];
	if (g->type != leaf)
		children[nc++] = g->child[0];
	if (g->type == divergence)
		children[nc++] = g->child[1];

	for (int i = 0; i < nc; i++)
		if (depth[children[i]->id] == -1)
		{
			parent[children[i]->id] = g;
			decompose_bc_helper(genes, ngenes, is_rbc, depth, low, parent, children[i], d+1);
			//Note: by default the DFS search starts at gene id 0 which is conventionally a leaf and therefore NOT an articulation point
			/*if (low[children[i]->id] >= depth[g->id] && depth[g->id] > 0)
			{
				//Is articulation node
				cout << g->id << " is an articulation node" << endl;
			}*/
			if (low[children[i]->id] > depth[g->id])
			{
				//(children[i], g) is an articulation edge
				//cout << "(" << g->id << "," << children[i]->id << ") is articulation edge" << endl;

				if (g->type != leaf && children[i] == g->child[0])
					is_rbc[children[i]->id] = true;
				else if (g->type == divergence && children[i] == g->child[1])
					is_rbc[children[i]->id] = true;
				else
					is_rbc[g->id] = true;
			}

			low[g->id] = min(low[g->id], low[children[i]->id]);
		}
		else if (children[i] != parent[g->id] || seen_parent)
			low[g->id] = min(low[g->id], depth[children[i]->id]);
		else		//Double edges between parent and child
			seen_parent = true;
}

void calculate_mpnsr(node** genes, int ngenes, node** species, int nspecies, node** species_map, reconciliation &r)
//Calculate MPNSR according to flow algorithm
//Assumes species tree is in bottom-up order
//Not complete (or even started)
{
	calculate_lca(genes, ngenes, species_map, r);

	//Convert to no-S reconciliation
	for (int g = 0; g < ngenes; g++)
		if (r.e[g] = S)
			r.e[g] = D;

	for (int s = 0; s < nspecies; s++)
	{
		//Determine set to move up
		//Move them up
		//r.f[i] = r.f[i]->parent[0];
	}
}

struct map
{
	set<edge> E;		//Note: using set<edge> instead of custom-built class assumes no double edges in gene network!
	int cost;
	map* back[2];
} ;

void merge(const set<edge> &e1, const set<edge> &e2, set<edge> &result)
{
	result = e1;

	edge sibling;

	for (set<edge>::iterator it = e2.begin(); it != e2.end(); it++)
	{
		if (it->child->type != root && it->parent->type == divergence)
		{
			sibling.parent = it->parent;

			if (it->parent->child[0] == it->child)
				sibling.child = it->parent->child[1];
			else
				sibling.child = it->parent->child[0];

			if (e1.find(sibling) != e1.end())
			{
				result.insert(edge(it->parent->parent[0], it->parent));
				result.erase(sibling);
				continue;
			}
		}

		result.insert(*it);
	}
}

ostream& operator<<(ostream &strm, const map &m)
{
	for (set<edge>::iterator it = m.E.begin(); it != m.E.end(); it++)
		strm << (*it) << ",";
	strm << " cost " << m.cost;
	return strm;
}

void calculate_mpr(node** genes, int ngenes, node** species, int nspecies, node** species_map, reconciliation &r)
//Calculate MPR according to new algorithm
{
	reconciliation lca, lca_hca;
	node* rbc[MAXID];
	int bc_num[MAXID];
	bool is_treechild[MAXID];
	int nbc;	

	//For now the max limit for each reticulation node is the LCA of the root of its biconnected component
	decompose_bc(genes, ngenes, rbc, bc_num, is_treechild, nbc);
	calculate_lca(genes, ngenes, species_map, lca);
	calculate_lca_hca(genes, ngenes, species_map, lca_hca);
	
	//Costs of reconciliations of BCC roots - we know they are mapped to LCAs, calculated above
	int bcc_root_cost[MAXID];

	//Initialise with leaves
	for (int i = 0; i < ngenes; i++)
	{
		if (genes[i]->type == leaf)
		{
			r.f[i] = species_map[i];
			r.e[i] = C;
			bcc_root_cost[i] = 0;
		}
		else
		{
			r.f[i] = NULL;
			bcc_root_cost[i] = INT_MAX;
		}
	}

	map active[MAXID][MAXID];		//May need to change second MAXID to higher
		//First index is species id, second index is edge set
	int nactive[MAXID];		//Index corresponds to first index of <active> (species id)
	node* root_bcc;
	int bcc_pointer = 0;

	//For each non-leaf BCC
	for (int bcc = 0; bcc < nbc; bcc++)
	{
		//If BCC is tree-child
		if (is_treechild[bcc])
		{
			int cost = 0, change;

			for (int i = 0; i < ngenes; i++)
				if (bc_num[genes[i]->id] == bcc)
				{
					r.f[i] = lca_hca.f[i];
					r.e[i] = lca_hca.e[i];

					//Calculate and add cost of edges to children
					cost += dist(r.f[i], r.f[genes[i]->child[0]->id]) * Lcost;
					if (r.e[i] == S)
						cost -= Lcost;
					if (bc_num[genes[i]->child[0]->id] != bcc)
						cost += bcc_root_cost[genes[i]->child[0]->id];

					if (genes[i]->type == divergence)
					{
						cost += dist(r.f[i], r.f[genes[i]->child[1]->id]) * Lcost;
						if (r.e[i] == S)
							cost -= Lcost;
						if (bc_num[genes[i]->child[1]->id] != bcc)
							cost += bcc_root_cost[genes[i]->child[1]->id];
					}

					//Add cost of duplication
					if (r.e[i] == D)
						cost += Dcost;
					
					//If root of BCC, then set cost and stop
					if (genes[i] == rbc[genes[i]->id])
					{
						bcc_root_cost[i] = cost;
						cout << "cost of BCC with root " << genes[i]->id << " is " << cost << endl;
						break;
					}
				}

			continue;
		}

		for (; bcc_pointer < ngenes; bcc_pointer++)		//Find root of current BCC - hacky
			if (bc_num[genes[bcc_pointer]->id] == bcc)
			{
				root_bcc = rbc[genes[bcc_pointer]->id];
				cout << "root of BCC " << bcc << " is " << root_bcc->id << endl;
				break;
			}

		//If BCC is not tree-child
		for (int sp = 0; sp < nspecies; sp++)
		{
			cout << "species " << sp << endl;

			nactive[sp] = 0;

			set<edge> Ex;

			//Find out-edges from current BCC that are mapped to current sp
			for (int i = 0; i < ngenes; i++)
			{
				if (genes[i] == root_bcc)
					break;

				//cout << "looking at " << i << endl;
				//cout << genes[i]->id << '\t' << rbc[genes[i]->id]->id << '\t' << bc_num[genes[i]->parent[0]->id] << '\t' << bcc << '\t' << lca.f[genes[i]->id]->id << '\t' << species[sp]->id << endl;

				if (genes[i] == rbc[genes[i]->id] && bc_num[genes[i]->parent[0]->id] == bcc && lca.f[genes[i]->id] == species[sp])	
					//Is root of BCC and parent is current BCC and LCA is current species
					Ex.insert(edge(genes[i]->parent[0], genes[i]));
			}

			//Initialise
			if (species[sp]->type == leaf)
			{
				active[sp][nactive[sp]].E = Ex;
				active[sp][nactive[sp]].cost = 0;
				for (set<edge>::iterator it = Ex.begin(); it != Ex.end(); it++)
					active[sp][nactive[sp]].cost += bcc_root_cost[it->child->id];
				active[sp][nactive[sp]].back[0] = active[sp][nactive[sp]].back[1] = NULL;
				nactive[sp]++;
			}
			//Or merge
			else
			{
				//Merge active[species[sp]->child[0]->id][], active[species[sp]->child[1]->id][], Ex
				int sp_l = species[sp]->child[0]->id, sp_r = species[sp]->child[1]->id;

				for (int i = 0; i < nactive[sp_l]; i++)
					for (int j = 0; j < nactive[sp_r]; j++)
					{
						//Merge feeders from inside BCC
						merge(active[sp_l][i].E, active[sp_r][j].E, active[sp][nactive[sp]].E);

						active[sp][nactive[sp]].cost = active[sp_l][i].cost + active[sp_r][j].cost;
						active[sp][nactive[sp]].cost += ( 2*active[sp][nactive[sp]].E.size() - active[sp_l][i].E.size() - active[sp_r][j].E.size() ) * Lcost;

						//Add existing out-edges
						active[sp][nactive[sp]].E.insert(Ex.begin(), Ex.end());

						for (set<edge>::iterator it = Ex.begin(); it != Ex.end(); it++)
							active[sp][nactive[sp]].cost += bcc_root_cost[it->child->id];

						//Set back-tracking
						active[sp][nactive[sp]].back[0] = &active[sp_l][i];
						active[sp][nactive[sp]].back[1] = &active[sp_r][j];

						nactive[sp]++;
					}
			}

			cout << "before expansion" << endl;
			for (int i = 0; i < nactive[sp]; i++)
				cout << active[sp][i] << " | ";
			cout << endl;

			edge sibling;
			map new_map;
			bool changed;
			
			//Expand reticulations and contract duplications
			for (int i = 0; i < nactive[sp]; i++)
			{
				changed = false;

				//Find smallest set and put in i - super slow and stupid but it should work
				int min_set = i;
				for (int j = i; j < nactive[sp]; j++)
					if (active[sp][j].E < active[sp][min_set].E)
						min_set = j;

				if (min_set != i)
				{
					//Swap min_set and i
					active[sp][nactive[sp]].E = active[sp][i].E;
					active[sp][i].E = active[sp][min_set].E;
					active[sp][min_set].E = active[sp][nactive[sp]].E;

					active[sp][nactive[sp]].cost = active[sp][i].cost;
					active[sp][i].cost = active[sp][min_set].cost;
					active[sp][min_set].cost = active[sp][nactive[sp]].cost;

					active[sp][nactive[sp]].back[0] = active[sp][i].back[0];
					active[sp][i].back[0] = active[sp][min_set].back[0];
					active[sp][min_set].back[0] = active[sp][nactive[sp]].back[0];

					active[sp][nactive[sp]].back[1] = active[sp][i].back[1];
					active[sp][i].back[1] = active[sp][min_set].back[1];
					active[sp][min_set].back[1] = active[sp][nactive[sp]].back[1];
				}

				//Finish if set contains root - since it is the smallest, there is no other set
				//Should it be root of BCC??
				if (active[sp][i].E.size() > 0 && active[sp][i].E.begin()->child == root_bcc)
					break;

				//cout << "considering: " << active[sp][i] << endl;

				//Contract duplications and expand forced reticulations in active[sp][i]
				//Remember there is no ordering in active[sp]
				for (set<edge>::iterator it = active[sp][i].E.begin(); it != active[sp][i].E.end(); it++)
				{
					if (it->parent->type == divergence)
					{
						sibling.parent = it->parent;

						if (it->parent->child[0] == it->child)
							sibling.child = it->parent->child[1];
						else
							sibling.child = it->parent->child[0];

						set<edge>::iterator sib_it = active[sp][i].E.find(sibling);

						if (sib_it != active[sp][i].E.end())
						{
							/*cout << "duplication" << endl;
							cout << active[sp][i] << endl;*/

							active[sp][i].E.insert(edge(it->parent->parent[0], it->parent));
							active[sp][i].E.erase(sib_it);
							active[sp][i].E.erase(it);

							active[sp][i].cost += Dcost;

							//Back-tracking - nothing changes

							changed = true;

							//cout << active[sp][i] << endl;
						}
					}
					else if (it->parent->type == reticulation)
					{
						//Forced reticulation
						if (lca.f[rbc[it->parent->id]->id] == species[sp])		//Change to better bound later
						{
							cout << "forced reticulation " << (*it) << endl;
							cout << active[sp][i] << endl;

							active[sp][i].E.insert(edge(it->parent->parent[0], it->parent));
							active[sp][i].E.insert(edge(it->parent->parent[1], it->parent));
							active[sp][i].E.erase(it);

							//Back-tracking - again nothing changes

							changed = true;

							cout << active[sp][i] << endl;
						}
					}
					
					//Only search ahead for matching sets
					if (changed)
					{

						for (int j = i+1; j < nactive[sp]; j++)
							if (active[sp][j].E == active[sp][i].E)
							{
								//Change back-tracking if j is cheaper than i, otherwise don't
								if (active[sp][j].cost < active[sp][i].cost)
								{
									active[sp][i].back[0] = active[sp][j].back[0];
									active[sp][i].back[1] = active[sp][j].back[1];
								}

								//Merge j and i into i
								active[sp][i].cost = min(active[sp][i].cost, active[sp][j].cost);
								nactive[sp]--;

								//Erase active[sp][j], replace with last element of active[sp]
								active[sp][j].E = active[sp][nactive[sp]].E;
								active[sp][j].cost = active[sp][nactive[sp]].cost;
								active[sp][j].back[0] = active[sp][nactive[sp]].back[0];
								active[sp][j].back[1] = active[sp][nactive[sp]].back[1];

								break;
							}

						break;
					}

				}

				//If changed, process again
				if (changed)
				{
					i--;
					continue;
				}

				//Expand non-forced reticulations (all, but one at a time) in active[sp][i]
				//But do not overwrite active[sp][i], put at the back of active[sp]
				for (set<edge>::iterator it = active[sp][i].E.begin(); it != active[sp][i].E.end(); it++)
				{
					if (it->parent->type == reticulation)
					{
						/*cout << "non-forced reticulation" << endl;
						cout << active[sp][i] << endl;*/

						active[sp][nactive[sp]].E = active[sp][i].E;

						active[sp][nactive[sp]].E.insert(edge(it->parent->parent[0], it->parent));
						active[sp][nactive[sp]].E.insert(edge(it->parent->parent[1], it->parent));
						active[sp][nactive[sp]].E.erase((*it));

						active[sp][nactive[sp]].cost = active[sp][i].cost;

						//Back-tracking
						active[sp][nactive[sp]].back[0] = active[sp][i].back[0];
						active[sp][nactive[sp]].back[1] = active[sp][i].back[1];

						//cout << active[sp][nactive[sp]] << endl;

						//Only search ahead
						for (int j = i+1; j < nactive[sp]; j++)
							if (active[sp][j].E == active[sp][nactive[sp]].E)
							{
								//If i is cheaper, change back-tracking, otherwise don't
								if (active[sp][i].cost < active[sp][j].cost)
								{
									active[sp][j].back[0] = active[sp][i].back[0];
									active[sp][j].back[1] = active[sp][i].back[1];
								}

								//Merge j and nactive[sp] into j
								active[sp][j].cost = min(active[sp][i].cost, active[sp][j].cost);
								nactive[sp]--;
								break;
							}

						nactive[sp]++;
					}
				}
			}

			cout << "after expansion" << endl;
			for (int i = 0; i < nactive[sp]; i++)
				cout << active[sp][i] << " | ";
			cout << endl;

			//At end, need to set bcc_root_cost
			//Triggers if current species is the LCA of the root of the BCC
			if (species[sp] == lca.f[root_bcc->id])
			{
				bcc_root_cost[root_bcc->id] = active[sp][0].cost;		//Assuming it has condensed. Should check

				//Resolve back-tracking
				resolve_backtracking(active[sp][0], species[sp], bcc, bc_num, r);

				//Not sure if we need this - maybe only for size 1 bcc
				//Also need to set r.e
				r.f[root_bcc->id] = lca.f[root_bcc->id];
				
				break;
			}
		}
	}

	cout << "MPR cost is " << bcc_root_cost[ngenes-1] << endl << endl;
}

void resolve_backtracking(map &top, node *species, int bcc, int* bc_num, reconciliation &r)
{
	cout << "backtracking " << top << " at " << species->id << endl;
	for (set<edge>::iterator it = top.E.begin(); it != top.E.end(); it++)
		if (bc_num[it->child->id] == bcc)
			set_from(it->child, top.back, species, bc_num, r);

	if (species->type == leaf)
		return;

	//Recurse
	if (top.back[0]->E.size() > 0)
		resolve_backtracking(*top.back[0], species->child[0], bcc, bc_num, r);
	if (top.back[1]->E.size() > 0)
		resolve_backtracking(*top.back[1], species->child[1], bcc, bc_num, r);
}

void set_from(node *gene, map **back, node *species, int *bc_num, reconciliation &r)
{
	if (species->type == leaf)
	{
		cout << "set " << gene->id << " at " << species->id << endl;

		r.f[gene->id] = species;

		switch (gene->type)
		{
			case reticulation: 
				r.e[gene->id] = R;
				if (r.f[gene->child[0]->id] == NULL && bc_num[gene->child[0]->id] == bc_num[gene->id])
					set_from(gene->child[0], back, species, bc_num, r);
				break;

			case root:
				r.e[gene->id] = T;
				if (r.f[gene->child[0]->id] == NULL && bc_num[gene->child[0]->id] == bc_num[gene->id])
					set_from(gene->child[0], back, species, bc_num, r);
				break;

			case divergence:
				r.e[gene->id] = D;
				if (r.f[gene->child[0]->id] == NULL && bc_num[gene->child[0]->id] == bc_num[gene->id])
					set_from(gene->child[0], back, species, bc_num, r);
				if (r.f[gene->child[1]->id] == NULL && bc_num[gene->child[1]->id] == bc_num[gene->id])
					set_from(gene->child[1], back, species, bc_num, r);
				break;
		}
	}
	else
	{
		//If node should be set in lower edge
		if (back[0]->E.find(edge(gene->parent[0], gene)) != back[0]->E.end()
				|| ( gene->type == reticulation && back[0]->E.find(edge(gene->parent[1], gene)) != back[0]->E.end() )
				|| back[1]->E.find(edge(gene->parent[0], gene)) != back[1]->E.end()
				|| ( gene->type == reticulation && back[1]->E.find(edge(gene->parent[1], gene)) != back[1]->E.end() ) )
			return;

		cout << "set " << gene->id << " at " << species->id << " above " << *back[0] << " | " << *back[1] << endl;

		r.f[gene->id] = species;

		switch(gene->type)
		{
			case reticulation:
				r.e[gene->id] = R;
				if (back[0]->E.find(edge(gene, gene->child[0])) == back[0]->E.end()
						&& back[1]->E.find(edge(gene, gene->child[0])) == back[1]->E.end()
						&& r.f[gene->child[0]->id] == NULL
						&& bc_num[gene->child[0]->id] == bc_num[gene->id])
					set_from(gene->child[0], back, species, bc_num, r);
				break;
			
			case root:
				r.e[gene->id] = T;
				//No need to recurse because root is always its own BCC
				break;

			case divergence:
				if (back[0]->E.find(edge(gene, gene->child[0])) != back[0]->E.end()
						&& back[1]->E.find(edge(gene, gene->child[1])) != back[1]->E.end())
				{
					r.e[gene->id] = S;
						return;
				}
				else if (back[1]->E.find(edge(gene, gene->child[0])) != back[1]->E.end()
						&& back[0]->E.find(edge(gene, gene->child[1])) != back[0]->E.end())
				{
					r.e[gene->id] = S;
					return;
				}
				else
				{
					r.e[gene->id] = D;

					if (r.f[gene->child[0]->id] == NULL
						&& bc_num[gene->child[0]->id] == bc_num[gene->id])
						set_from(gene->child[0], back, species, bc_num, r);

					if (r.f[gene->child[1]->id] == NULL
						&& bc_num[gene->child[1]->id] == bc_num[gene->id])
						set_from(gene->child[1], back, species, bc_num, r);
				}

				//Weird logic but should work
				//Need re-assessing
				/*if (back[0]->E.find(edge(gene, gene->child[1])) == back[0]->E.end()
						&& back[1]->E.find(edge(gene, gene->child[1])) == back[1]->E.end()
						&& r.f[gene->child[1]->id] == NULL
						&& bc_num[gene->child[1]->id] == bc_num[gene->id])
					set_from(gene->child[1], back, species, bc_num, r);*/
				break;
		}
	}
}

void create_edge_list(node** genes, int ngenes, edge** edges, int &nedges)
	//Again unused and needs changing
{
	for (int i = 0; i < ngenes; i++)
	{
		edges[nedges] = new edge;
		edges[nedges]->child = genes[i];
		edges[nedges]->parent = (genes[i]->type == root ? NULL : genes[i]->parent[0]);
		nedges++;

		if (genes[i]->type == reticulation)
		{
			edges[nedges] = new edge;
			edges[nedges]->child = genes[i];
			edges[nedges]->parent = genes[i]->parent[1];
			nedges++;
		}
	}
}

/*To do (28/8/19):
 * Fix iteration of active set		./
 * File input for gene/species network/trees		./
 * More testing					./
 * LCA-HCA for tree-child BCCs		./
 * Back-tracking to generate MPR		./
 * MPNSR algorithm
 * Simulation to generate gene/species network/trees (separate program)
 */
