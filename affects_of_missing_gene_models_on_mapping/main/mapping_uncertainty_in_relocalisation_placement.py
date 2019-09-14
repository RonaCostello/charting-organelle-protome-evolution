import pandas as pd

# def main():
#     with open('../data/Phyldog_relocalisation_duplication_data.csv') as input_file:
#         for line in input_file:
#             token = line.rstrip().split()
#             print(token)
#
# if __name__ == '__main__':
# 	main()
#
#
#
#
#QUERY N3 AS TEST TO SEE HOW DIFF TO ORIGINAL SCRIPTTT!!!!!!!!!!!!!!

# query the relocalisation dataframe:
df = pd.read_csv('../data/Phyldog_relocalisation_duplication_data.csv')

df['C_ingroup_retention'] = pd.to_numeric(df["C_ingroup_retention"], errors='coerce')
df['C_outgroup_retention'] = pd.to_numeric(df["C_outgroup_retention"], errors='coerce')

print(df.loc[(df['species_tree_node'] == 'N3') & (df['C_relocalisation'] == 'gain') & (df['C_ingroup_retention'] > 0.75) & (df['C_outgroup_retention'] >= 0.75)])
