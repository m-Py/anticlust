int: number_nodes;
int: number_clusters;
int: number_edges;

set of int: iteration_nodes = 1..number_nodes;
set of int: iteration_clusters = 1..number_clusters;
set of int: iteration_edges = 1..number_edges;

array[iteration_clusters] of int: target_groups;
array[iteration_edges] of int: edges1;
array[iteration_edges] of int: edges2;

%array[1..number_nodes] of var iteration_clusters: colors_of_nodes;
array[iteration_nodes] of var int: colors_of_nodes;
constraint forall(i in iteration_nodes) (colors_of_nodes[i] in iteration_clusters);

constraint forall(i in iteration_edges)(colors_of_nodes[edges1[i]] != colors_of_nodes[edges2[i]]);
constraint forall(i in iteration_clusters)(count(colors_of_nodes, i) <= target_groups[i]);
%constraint symmetry_breaking_constraint(colors_of_nodes[1] == 1);