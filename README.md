# coexpression_convergence

## 1. You can access the CommonMind Consortium RNA sequencing data from Synapse. 

You can also use other transcriptmic datasets (i.e. GTEx). 
The data will need to formatted such that each column represents a gene and each row represents an individual. 



## 2. Calculate coexpression using Calculate.coexpression.R

This will calculate the Pearson's pairwise coexpression data.



## 3. To run convergence permutation, use Permutation.R 

Run the following command: 

```
Rscript Permutation.R <genelist.txt> <# of permutations> <phenoname> 
```

genelist.txt is a list of ensembl IDs for your genes of interest. 

The # of permutations is how many permutations you want to run. We would suggest at least 1,000,000 permutations. 

phenoname is what you want to name your output


