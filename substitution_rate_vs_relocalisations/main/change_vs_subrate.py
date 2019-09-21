from ete3 import Tree
import pandas as pd
import csv
import itertools as it

def branch_length_between_nodes(tree, node, parent_node):
    start_node = tree&node
    end_node = tree&parent_node
    common_ancestor = (tree.get_common_ancestor(start_node, end_node))
    if common_ancestor != end_node and start_node != end_node:
        raise Exception('node2 must be an ancestor of node1 and not a sister node')

    return(tree.get_distance(node, parent_node))

def is_relocalisation_branch(row, df):
    evidence_for_relocalisation = False
    high_scoring_relocalisation = False

    for i, col in enumerate(df.columns):
        if getattr(row, col) == 'gain' or getattr(row, col) == 'loss':
            evidence_for_relocalisation = True
            ingroup_retention = float(getattr(row, (df.columns)[i+1]))
            outgroup_retention = float(getattr(row, (df.columns)[i+2]))
            if ingroup_retention >= 0.75 and outgroup_retention >= 0.75:
                high_scoring_relocalisation = True

    return(evidence_for_relocalisation, high_scoring_relocalisation)

def main():
    results_relocalisation_branches = []
    nodes_relocalisation_branches = []
    results_nonrelocalisation_branches = []
    nodes_nonrelocalisation_branches = []

    species_tree = Tree("../data/Phytozome10_constrainedTree_rooted_labelled.tree", format=1)

    input_df = pd.read_csv('../data/Phyldog_relocalisation_duplication_data.csv')
    cols=[i for i in input_df.columns if 'retention' in i]
    for col in cols:
        input_df[col]=pd.to_numeric(input_df[col], errors='coerce')

    speciation_df = input_df.loc[(input_df['orthogroup_tree_node_event'] == 'Speciation') & (input_df['retention_score'] == 100.0)]
    for row in speciation_df.itertuples(index=False):
        orthogroup = getattr(row, "orthogroup")
        gene_tree = Tree(f'../data/Phyldog_trees/{orthogroup}.locus.tree', format=1)

        gene_tree_node = getattr(row, "orthogroup_tree_node")
        gene_tree_node_mapped = getattr(row, "species_tree_node")

        gene_tree_node_parent = (gene_tree&gene_tree_node).up.name

        if ((gene_tree&gene_tree_node).up).is_root():
            continue
        else:
            gene_tree_node_parent_mapped = input_df.loc[(input_df['orthogroup'] == orthogroup) & (input_df['orthogroup_tree_node'] == gene_tree_node_parent), 'species_tree_node'].iloc[0]

        if ((species_tree&gene_tree_node_mapped).up) == species_tree&gene_tree_node_parent_mapped:
            species_tree_distance = float(branch_length_between_nodes(species_tree, gene_tree_node_mapped, gene_tree_node_parent_mapped))
            gene_tree_distance = float(branch_length_between_nodes(gene_tree, gene_tree_node, gene_tree_node_parent))
            substitution_rate = gene_tree_distance/species_tree_distance
            is_branch_reloc, is_branch_retained_reloc = is_relocalisation_branch(row, input_df)
            if is_branch_reloc == True:
                if is_branch_retained_reloc == True:
                    results_relocalisation_branches.append(substitution_rate)
                    nodes_relocalisation_branches.append(gene_tree_node_mapped)
                else:
                    continue
            else:
                results_nonrelocalisation_branches.append(substitution_rate)
                nodes_nonrelocalisation_branches.append(gene_tree_node_mapped)

    with open('substitution_rate_on_relocalisation_branches_with_nodes.csv', 'w') as f:
        csv.writer(f).writerow(["substitution_rate_on_relocalisation_branches", "species_tree_node_for_reloc_branch" "substitution_rate_on_non_relocalisation_branches", "species_tree_node_for_nonreloc_branch"])
        csv.writer(f).writerows(it.zip_longest(results_relocalisation_branches, nodes_relocalisation_branches, results_nonrelocalisation_branches, nodes_nonrelocalisation_branches))



if __name__ == '__main__':
	main()
