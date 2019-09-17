from ete3 import Tree

def topology_between_nodes(tree, node1, node2): # node1 is the node closest to the terminal brances of the tree, with node2 being the deeper node in the tree
    start_node = tree&node1
    end_node = tree&node2

    distance = tree.get_distance(start_node, end_node, topology_only=True)

    nodes_between = []
    for i in range(int(distance)):
        nodes_between.append((start_node.up).name)
        start_node = start_node.up

    return(nodes_between)

def branch_length_between_nodes(tree, node1, node2):
    all_nodes = topology_between_nodes(tree, node1, node2)
    all_nodes.append(node1)
    print(all_nodes)

def main():
    species_tree = Tree("../data/Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
    branch_length_between_nodes(species_tree, 'N18', 'N18')

if __name__ == '__main__':
	main()
