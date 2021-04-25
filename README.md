# pancreatic_cancers
The objective of this repository is to explore working with gene expression datasets, taking the use case of Pancreatic Adenocarcinoma.

## Dataset
1. We have a dataset containing the fene expression data for pancreatinc cancer samples. The data consists of 18465 genes and 183 samples.
2. 25 genes list defining the IFN signature
3. Gene expression in a Normal Pancrease tissue

### Task 1 - Cleaning the dataset, generating gene expression distribution for all samples
1. Removing genes with Null values - **4367** genes have Null values
2. Of the remaining genes, we sort the dataframe based on maximum gene sum values and create a heatmap (using the seaborn package)
![](https://github.com/raunak95/pancreatic_cancers/blob/master/results/task1/gene_distribution.png)

### Task 2 - Subsetting the data for only Exocrine tumors and removing Neuroendocrine tumors using PCA.
1. Using the metadata information to separate the two types of tumors, and running PCA on the data using only the first two Prinicipal components - the two types of tumors are fairly separable
![](https://github.com/raunak95/pancreatic_cancers/blob/master/results/task2/PCA_scatter_tumor_type.png)
