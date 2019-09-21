import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np

data = pd.read_csv('substitution_rate_on_relocalisation_branches.csv')

relocalisation_branch_sub_rate = data.relocalisation_branches.tolist()
cleaned_relocalisation_list = [x for x in relocalisation_branch_sub_rate if str(x) != 'nan']
non_relocalisation_branch_sub_rate = data.non_relocalisation_branches.tolist()

relocalisation_log10 = [math.log10(number) for number in cleaned_relocalisation_list]
non_relocalisation_log10 = [math.log10(number) for number in non_relocalisation_branch_sub_rate]

combined_list = []
combined_list.append(relocalisation_log10)
combined_list.append(non_relocalisation_log10)

# plt.boxplot(combined_list)
# plt.savefig('foo2.pdf')



plt.hist(relocalisation_log10, bins=30)
plt.ylabel('Relocalisation_branches')
plt.savefig('hist_relocs.pdf')

plt.hist(non_relocalisation_log10, bins=30)
plt.ylabel('Non_relocalisation_branches')
plt.savefig('hist_Non_relocs.pdf')


# x = [[1.2, 2.3, 3.0, 4.5],
#      [1.1, 2.2, 2.9]]
# plt.boxplot(x)
# plt.show()
