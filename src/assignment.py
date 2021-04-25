import os, time
from cmapPy.pandasGEXpress import parse_gct
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd

main_dir = '/home/raunak/Projects/Bioinformatics/ElucidataAssignment'

def sortDF(df):
    df2 = df.T
    df2['total'] = df2.sum(axis=1)
    df2 = df2.sort_values(['total'], ascending = False).drop(['total'], axis=1).T
    return df2

data = parse_gct.parse(os.path.join(main_dir, 'data/cancer_dataset.gct'))
df = data.data_df.T
df_meta_col = data.col_metadata_df
df_meta_row = data.row_metadata_df

## Task 1
# Remove genes with NANs
df2 = df.dropna(axis = 1)
df2 = sortDF(df2)
missing_genes = list(set(df.columns) - set(df2.columns))
print('Number of genes having NaNs = {}'.format(len(missing_genes)))

print('Generating heatmap for gene distribution')
sns.heatmap(df2)
plt.gcf().set_size_inches((25, 25))
plt.savefig(os.path.join(main_dir, 'results/task1/gene_distribution.png'))
plt.close('all')

## Task 2
print('Running PCA')
pca = PCA(n_components=2)
pca_data = pca.fit_transform(df2)
df_pca = pd.DataFrame(data = pca_data, columns = ['PC1', 'PC2'], index = df2.index)
df_pca.to_csv(os.path.join(main_dir, 'results/task2/PCA_output.csv'))
plt.scatter(df_pca['PC1'], df_pca['PC2'])
plt.gcf().set_size_inches((25, 25))
plt.savefig(os.path.join(main_dir, 'results/task2/PCA_scatter.png'))
plt.close('all')

# Overlay neuroendocrine tumour info
df_meta_pca = df_meta_col[['histological_type_other']].fillna('')
df_pca2 = pd.merge(df_pca, df_meta_pca, left_index = True, right_index = True)
df_pca2['tumor_type'] = df_pca2['histological_type_other'].map(lambda x: 'Neuroendocrine' if 'neuroendocrine' in x.lower() else 'Exocrine')

groups = df_pca2.groupby('tumor_type')

for name, group in groups:
    plt.plot(group['PC1'], group['PC2'], marker="o", linestyle="", label=name)

plt.legend()
plt.gcf().set_size_inches((10, 10))
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PC1 v/s PC2 for Pancreatic tumor types')
plt.savefig(os.path.join(main_dir, 'results/task2/PCA_scatter_tumor_type.png'))
plt.close('all')

# Number of components vs captured variance
print('Calculating no. of components v/s captured variance')
captured_variance = {}
st = time.time()
for i in range(2, 51):
    pca = PCA(n_components=i)
    pca.fit_transform(df2)
    captured_variance[i] = pca.explained_variance_ratio_.cumsum()[-1]
print('Time taken = {} seconds'.format(time.time() - st))

for n in [2, 10, 20, 30, 50]:
    print('When PCA components = {}, the captured variance = {}'.format(n, captured_variance[n]))

df_pca_var = pd.DataFrame.from_dict(captured_variance, orient = 'index')
df_pca_var.plot(legend = None)
plt.gcf().set_size_inches((10, 10))
plt.xlabel('No. of components')
plt.ylabel('Captured variance')
plt.title('Captured variance versus no. of PCA components')
plt.savefig(os.path.join(main_dir, 'results/task2/variance_nComponents.png'))
plt.close('all')

## Task 3
# Separate out exocrine tumors
df_exocrine = df_pca2[df_pca2['tumor_type'] == 'Exocrine']
ids_exocrine = df_exocrine.index.tolist()
df_exo = df2[df2.index.isin(ids_exocrine)].T
# Save df_exo and run GSVA
'''
Command on Unix:
docker run -v $(pwd):$(pwd) vacation/gsva:1.0.4 GSVA --gmt $(pwd)/gene_set.gmt $(pwd)/exocrine_mat.csv --output $(pwd)/pathways_exocrine.csv
'''
df_path_exo = pd.read_csv(os.path.join(main_dir, 'data/pathways_exocrine.csv')).drop(['name'], axis=1).T
df_path_exo.columns = ['score']
plt.hist(df_path_exo['score'].tolist(), bins=40)
plt.gcf().set_size_inches((15, 15))
plt.title('GSVA scores historgram for Exocrine tumors')
plt.savefig(os.path.join(main_dir, 'results/task3/GSVA_historgram.png'))
plt.close('all')

## Bonus task
# GSVA analysis for pancreatic_cancers
df_pan = pd.read_csv(os.path.join(main_dir, 'results/bonus/pancreatic_cancer_pathways.csv')).drop(['name'], axis=1).T
df_pan.columns = ['score']
plt.hist(df_path_exo['score'].tolist(), bins=40)
plt.gcf().set_size_inches((15, 15))
plt.title('GSVA scores historgram for Pancreatic cancers')
plt.savefig(os.path.join(main_dir, 'results/bonus/GSVA_pancreaticCancer_historgram.png'))
plt.close('all')