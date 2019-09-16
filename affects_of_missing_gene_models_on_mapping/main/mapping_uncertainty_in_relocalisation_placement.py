import pandas as pd
from ete3 import Tree

def set_up_results_file(tree):
    organelles = ['Chloroplast', 'Mitochondria', 'Secretory', 'Peroxisome']
    results_df = pd.DataFrame(columns=['Species_tree_node', 'Organelle', 'old_gains', 'old_losses', 'certain_gains', 'certain_losses', 'new_gains', 'new_losses'])
    row = 0
    for organelle in organelles:
        for node in tree.traverse():
            results_df.loc[row] = [node.name, organelle] + [0]*6
            row+=1

    return(results_df)

def topology_between_nodes(tree, node1, node2):
    start_node = tree&node1
    end_node = tree&node2

    distance = tree.get_distance(start_node, end_node, topology_only=True)

    nodes_between = []
    for i in range(int(distance)):
        nodes_between.append((start_node.up).name)
        start_node = start_node.up

    return(float(distance), nodes_between)

def add_count_to_node(nodes, organelle, df, relocalisation, column_selection):
    columns_to_update = []
    cols=[i for i in df.columns if relocalisation in i]
    for col in cols:
        if (col.rstrip().rsplit('_'))[0] in column_selection:
            columns_to_update.append(col)

    increment = 1.0/len(nodes)
    print(increment)
    df.loc[(df['Species_tree_node'].isin(nodes)) & (df['Organelle'] == organelle), columns_to_update] += increment

    return(df)


def main():
    species_tree = Tree("../data/Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
    results_df = set_up_results_file(species_tree)

    input_df = pd.read_csv('../data/Phyldog_relocalisation_duplication_data.csv')

    cols=[i for i in input_df.columns if 'retention' in i]
    for col in cols:
        input_df[col]=pd.to_numeric(input_df[col], errors='coerce')
    for organelle in ['Chloroplast']:
        relocalisation_df = input_df.loc[(input_df[f'{organelle[0]}_ingroup_retention'] >= 0.75) & (input_df[f'{organelle[0]}_outgroup_retention'] >= 0.75)]

        for row in relocalisation_df.itertuples(index=False):
            orthogroup = getattr(row, "orthogroup")
            relocalisation_type = getattr(row, f"{organelle[0]}_relocalisation")
            gene_tree = Tree(f'../data/Phyldog_trees/{orthogroup}.locus.tree', format=1)

            gene_tree_node = getattr(row, "orthogroup_tree_node")
            species_tree_node = getattr(row, "species_tree_node")
            gene_tree_node_parent = (gene_tree&gene_tree_node).up.name
            species_tree_node_parent = input_df.loc[(input_df['orthogroup'] == orthogroup) & (input_df['orthogroup_tree_node'] == gene_tree_node), 'species_tree_node'].iloc[0]

            if 'N' not in species_tree_node_parent:
                results_df = add_count_to_node([species_tree_node_parent], organelle, results_df, relocalisation_type, ['old', 'certain', 'new'])

    results_df.to_csv('test.csv')
            # print(species_tree_node_parent)










    # distance, node_list = topology_between_nodes(species_tree, 'N35', 'N28')



if __name__ == '__main__':
	main()

#
#
#
#QUERY N3 AS TEST TO SEE HOW DIFF TO ORIGINAL SCRIPTTT!!!!!!!!!!!!!!

# query the relocalisation dataframe:
# df = pd.read_csv('../data/Phyldog_relocalisation_duplication_data.csv')
#
# df['C_ingroup_retention'] = pd.to_numeric(df["C_ingroup_retention"], errors='coerce')
# df['C_outgroup_retention'] = pd.to_numeric(df["C_outgroup_retention"], errors='coerce')
# print("loss")
# df = (df.loc[(df['species_tree_node'] == 'N13') & (df['C_relocalisation'] == 'loss') & (df['C_ingroup_retention'] >= 0.75) & (df['C_outgroup_retention'] >= 0.75)])
# df.to_csv("test.csv")
# print(df)
