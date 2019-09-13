<img src="https://i.imgur.com/sRIxQAf.png" alt="logo" width=200>

### Differential gene coexpression analysis and network model selection across multiple subgroups

SuperNOVA uses the robust ANOVA-based Chow test to determine whether two genes are better represented by a global correlation model or multiple correlation models between two or more subgroups of samples. 

<img src="https://i.imgur.com/vtRwYLN.png" alt="explain0" width=500>

### Benefits

**Compared to other existing tools, SuperNOVA can:**
* perform differential gene coexpression analysis on more than two subgroups at the same time, greatly increasing the power to detect differential correlation and simplifying downstream analysis 
* detect both changes in correlation strength and in correlation slope, identifying important changes in gene function and regulation that would otherwise be overlooked by existing tools
* generate subgroup-specific gene coexpression network models from the global model as a background, enabling robust subgroup coexpression network analyses even when the sample size is small. 
* evaluate trillions of features simultaneously using a parallel backend

SuperNOVA is the first differential correlation tool that can accurately identify changes in correlation slope, while being just as accurate as Fisher z-transformed differential correlation tools for changes in correlation strength.

<img src="https://i.imgur.com/ww9jlDt.png" alt="res" width=500>

### Install 

```R
devtools::install_github("ryananeff/superNOVA")
```

### Load superNOVA

```R
library("superNOVA")
```

### Prepare input expression data and design matrices

SuperNOVA needs an expression matrix and a design matrix to run.

Expression matrix: a features (rows) by samples (columns) matrix of continuous numerical values.
Design matrix: a samples (rows) by groups (columns) matrix, where a sample's membership in a group is represented by a binary numerical value of 1 (yes) or 0 (no).

NOTE: the samples in the expression matrix should exactly match the samples in the design matrix. 

```R
#load data
datExpr = read.delim(file = "expression_file.tsv",sep = "\t",header = T) #expression
design_mat = read.delim(file = "design_file.tsv", sep="\t", header=T) #design

# make model.matrix from design matrix, where group1 is the first group, group2 is the second group, etc...
design_mat = model.matrix(~0+group1+group2+group3+...,data=data.frame(design_mat))

```

### Run SuperNOVA

**Single threaded mode**

```R
# Run the analysis
supernova_obj = superNOVA::chowCor(design_mat=design_mat,
                              matA = datExpr,
                              compare = c("group1","group2","group3",...),
                              corrType="pearson")

# Flatten the intermediate object to a flat matrix
supernova_table = superNOVA::flattenChow(supernova_obj)

# Write output
write.table(supernova_table,file="superNOVA_results.tsv", sep = "\t",quote = F,row.names = F)

```

**Batch mode (using the [batchtools](https://mllg.github.io/batchtools/index.html) R package)**

```R
corrType = "pearson"
split = 10	                            #Number of times to split the input into smaller chunks
outputfile = "superNOVA_results.tsv"	  #output file (TSV)
chunkSize = 1                           #Increase this number if it takes a long time to submit jobs to HPC
batchDir = "supernova_batch_directory/" #where to store the intermediate results

walltime = 30 								#walltime in minutes
memory = 10000								#memory per core (MB)

chowParallel(datExpr, design_mat, outputfile, sigOutput=T, corrType=corrType,
             perBatch=split,timePerJob = walltime, memPerJob = memory,
             batchDir = batchDir, chunkSize=chunkSize)

```

### Output file structure

| Gene1                      | Gene2                       | groupCor                                                                        | groupCorPval                                                                                      | globalCor                                                    | globalCorP                                                                    | pValDiff                                                                                                                                                                                             | qValDiff                                                | Classes                                                                                  | group_order                                                          |
|----------------------------|-----------------------------|---------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|--------------------------------------------------------------|-------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|------------------------------------------------------------------------------------------|----------------------------------------------------------------------|
| The first gene in the pair | The second gene in the pair | The R^2 value between the pair for each subgroup of the data, separated by '/'. | The p-value of the correlation between the pair for each subgroup of the data, separated by '/'.  | The R^2 value between the pair for all of the data together. | The p-value of the correlation between the pair for all of the data together. | The p-value that a single linear regression model fits all of the data vs. separate subgroup models under ANOVA. | The adjusted q-value of the reported Chow test p-value. | The correlation direction ("+","0","-") for each subgroup of the data, separated by '/'. | The order of the subgroups evaluated by SuperNOVA, separated by '/'. |


### Where to find more help

More help can be found in the documentation for each function within R by typing in `??superNOVA`.


### 'I got an error message, what does this mean?'

After you have confirmed the inputs match the required inputs as per the documentation, submit a Github issue with your code, some example data that reproduces the issue, and the output from `sessionInfo()`, and I will look at it right away.

### Developed by 

Ryan Neff at the Bin Zhang lab, Icahn School of Medicine at Mount Sinai. 

Special thanks to: Andy McKenzie, Kevin Bu, Bin Zhang, Alice Wang, and everyone at the Zhang lab who has supported this project and given invaluable feedback.

### News

*9/13/2019:* Note that due to a bug in the R batchtools package currently, please either run this tool in a background thread (non-interactively) or, if you are running it interactively with screen, to detach the screen before it gets to the "Submitting..." progress bar and not reattach. If you are running it interactively outside of screen, do not resize the terminal window at all while it says "Submitting..." or "Waiting...". This is due to a bug in the progress bar which causes the batchtools package to crash unexpectedly as shown in the Github issue link above.

To mitigate this issue (which disables the progress bar completely), add this to the top of your R file: 
`options(batchtools.progress = FALSE)`
