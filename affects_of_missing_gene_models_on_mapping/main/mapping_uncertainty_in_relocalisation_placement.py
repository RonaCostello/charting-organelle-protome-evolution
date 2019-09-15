import pandas as pd
from ete3 import Tree

def topology_between_nodes(tree, node1, node2):
    start_node = tree&node1
    end_node = tree&node2

    distance = tree.get_distance(start_node, end_node, topology_only=True)

    nodes_between = []
    for i in range(int(distance)):
        nodes_between.append((start_node.up).name)
        start_node = start_node.up

    return(distance, nodes_between)

def main():
    species_tree = Tree("../data/Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
    


    distance, node_list = topology_between_nodes(species_tree, 'N35', 'N14')


    # with open('../data/Phyldog_relocalisation_duplication_data.csv') as input_file:
    #     for line in input_file:
    #         token = line.rstrip().split()
    #         print(token)

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
