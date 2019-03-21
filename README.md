
# Documentation and software for the "Predicting drug-pair synergy from the predicted synergy probabilities of individual drugs"

This README provides a high-level overview of the steps needed to reproduce the results in drug-pair synergy publication. The three primary steps are:

1. Download and process the synergy, gene expression and 2D drug structure data
2. Build the synergy prediction program
3. Run the synergy prediction program

Please let me know (jgans@lanl.gov) if you encounter any problems building or running these programs.

## Download and preprocess the data

### Cell line gene expression features

The cell line gene expression features are the RNA Seq RPKM values (provided by the CCLE) for the LINC1000 genes. The CCLE RNA Seq RPKM values are available [here](https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct); note that you will need to create an account and login to access the CCLE data. The data file that should be downloaded is called: `CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct`. (There is already a newer file of RNA Seq RPKM values in the CCLE data respository, so the file listed above is displayed under the "Previous Releases" header on the CCLE data webpage).

The perl script, `jdacs4c-staging/LANL_Synergy/preprocess/CCLE_format_rnaseq.pl`, is provided to format the CCLE data and extract the LINCS100 gene expression values corresponding to the LINCS1000 genes. After downloading the CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct file, run and redirect the script:

`CCLE_format_rnaseq.pl CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct > CCLE_DepMap_18q3_RNAseq_RPKM.LINCS1000.tsv`

to produce a tab-delimited file with cell-line samples corresponding to the rows and gene names corresponding to the column. Note that some of the LINCS1000 gene names are not found in the CCLE input file (which will generate `Failed to match LINCS gene` warnings. However, a total of 958 genes should be matched.

### 2D drug structure features

The drug features are the 1021-bit binary fingerprints computed from desalted, 2D chemical structres using the [OpenBabel program](http://openbabel.org/wiki/Main_Page). The 2D drug structures for both the NCI-ALMANAC and Merck data sets are contained in the `jdacs4c-staging/LANL_Synergy/preprocess/ALMANAC_and_Merck.smiles` SMILES file. The fingerprints, computed using the following commands: 

```
obabel ALMANAC_and_Merck.smiles -ofps -OALMANAC_and_Merck.obable.FP2 -xfFP2
cat ALMANAC_and_Merck.obable.FP2 | ./expand_fingerprints.pl > ALMANAC_and_Merck.obable.FP2.csv
```

are provided in the file `jdacs4c-staging/LANL_Synergy/preprocess/ALMANAC_and_Merck.obable.FP2.csv`. The `jdacs4c-staging/LANL_Synergy/preprocess/expand_fingerprints.pl` script converts the hex-encoded drug fingerprints into a text-based binary fingerprint (i.e. strings of 1's and 0's). Please note that CID drug IDs are used for both the NCI-ALMANAC and Merck drugs (the NSC drug IDs provide by the NCI are not used).

### NCI-ALMANAC drug synergy data
  
The NCI-ALMANAC drug-pair synergy data is available from [ComboDrugGrowth_Nov2017.zip](https://wiki.nci.nih.gov/download/attachments/338237347/ComboDrugGrowth_Nov2017.zip). Unzip this file to extract `ComboDrugGrowth_Nov2017.csv`. 

The first preprocessing step for the NCI-ALMANAC data is to (a) extract the subset of the data that has expression data in CCLE cell line gene expression data set and (b) converting the NSC-based drug ID used by NCI-ALMANAC to the CID-based drug ID used to label the drug features. Both of these steps are performed by the `jdacs4c-staging/LANL_Synergy/preprocess/ALMANAC_to_cid_and_CCLE.pl` perl script. This script contains the hard-coded mappings between NSC <-> CID drug ids and the CCLE <-> NCI60 cell line names. The script is run as:

`ALMANAC_to_cid_and_CCLE.pl ComboDrugGrowth_Nov2017.csv > ComboDrugGrowth_Nov2017.CID.CCLE.csv`

where the output (written to STDOUT) is redirected to a filename you select (I used `ComboDrugGrowth_Nov2017.CID.CCLE.csv` in the example above).

The next step of preprocessing the NCI-ALMANAC data is to compute the drug pair synergy values. Since the NCI-ALAMANC data file contains both the single agent and drug pair responses (for the different concentrations tested), parsing and processing this file is somewhat challenging. To manage the parsing complexity, this preprocessing step is performed by a C++ based program called `gemini_prep`. 

However, please note that the `gemini_prep` program requires the [GNU Scientific Library](https://www.gnu.org/software/gsl/), in addition to C++ compiler (the code has been tested with g++ ver. 4.4.7). After installing the GNU Scientific Library, you will need to edit `jdacs4c-staging/LANL_synergy/preprocess/Makefile` to specify the directories that contain the GSL library files (`GSL_LIB_DIR`) and GSL include files (`GSL_INCLUDE_DIR`). If your compiler does not support OpenMP (i.e. you have a Mac with the stock clang compiler as of 2018) then comment out the `-fopenmp` and the program should still compile (but run more slowly, since multi-threading is now disabled). 

The `gemini_prep` is compiled by running the `make` command in the `jdacs4c-staging/LANL_Synergy/preprocess` directory. 

After you run `make` (which will build both the `gemini_prep` and another program called `synergy_search` that will be described below), you will be able to run the `gemini_prep` command using:

```
./gemini_prep \
	-i /home/jgans/Cancer/data/ALMANAC/ComboDrugGrowth_Nov2017.CID.CCLE.csv \
	-p ALMANAC_CID_CCLE_
 ```

where the input file (passed as the argument to the `-i` flag) is the output of the `jdacs4c-staging/LANL_Synergy/preprocess/ALMANAC_to_cid_and_CCLE.pl` script and the argument to the `-p` flag specifies the output prefix that will be used for all of the output files (i.e. `ALMANAC_CID_CCLE_` in the above example). The `gemini_prep` program will generate a total of six output files (all starting with the specified prefix, i.e. `ALMANAC_CID_CCLE_`):

* `ALMANAC_CID_CCLE_bliss_ave_synergy.csv`
* `ALMANAC_CID_CCLE_bliss_stdev_synergy.csv`
* `ALMANAC_CID_CCLE_loewe_ave_synergy.csv`
* `ALMANAC_CID_CCLE_loewe_stdev_synergy.csv`
* `ALMANAC_CID_CCLE_pair_ave_growth.csv`
* `ALMANAC_CID_CCLE_pair_stdev_growth.csv`

The required synergy for training and testing machine learning models is contained in two of these files:

* `ALMANAC_CID_CCLE_bliss_ave_synergy.csv` for the Bliss synergy model
* `ALMANAC_CID_CCLE_loewe_ave_synergy.csv` for the Loewe synergy model

These files contain the average (over replicate experiments) minimum synergy (i.e. most synergistic over all tested concentrations).

### Merck drug synergy data

The Merck drug-pair synergy data is contained in the supplementary online data for O'Neil et. al. "An Unbiased Oncology Compound Screen to Identify Novel Combination Strategies", Molecular Cancer Therapeutics, 2016 Jun;15(6):1155-62. The single agent response data is stored in one [Excel file](http://mct.aacrjournals.org/highwire/filestream/53222/field_highwire_adjunct_files/1/156849_1_supp_0_w2lh45.xlsx) and the combination response data is stored in another [Excel file](http://mct.aacrjournals.org/highwire/filestream/53222/field_highwire_adjunct_files/3/156849_1_supp_1_w2lrww.xls).

To process the Merck data, both the single agent and combination response files must be manually converted to comma delimited files. The provided perl script, `jdacs4c-staging/LANL_Synergy/preprocess/process_merck.pl`, converts these CSV files (*not* the Excel files!) into a single file that contains the most synergistic measurements for each drug pair and cell line tested. This file is run as:

`./process_merck.pl --single <single agent response CSV file> --pair <combination response CSV file> > merck_pair.csv`

where the output (written to STDOUT) must be redirected to a filename you specify (i.e. `merck_pair.csv` in this example).

### Extract the per-cell line synergy values for both NCI-ALMANAC and Merck

The final preprocessing step for both the NCI-ALMANAC and the Merck data sets is to extract the synergy values for each cell line as a separate file. This task is performed by the `synergy_search` C++ program that should have been built when the `make` command was run in the `jdacs4c-staging/LANL_synergy/preprocess` directory (see the instructions above for the NCI-ALMANAC data). The process of creating the per-cell line drug-pair synergy files is automated by the provided shell script `batch_synergy_search.sh`. This shell script expects to create the output synergy files in two existing directories: `data_bliss/` and `data_loewe/`. If you would like to place these files in a different location, please edit the `OUTPUT_BLISS_DIR` and `OUTPUT_LOEWE_DIR` at the top of the script. This script is run as:

`./batch_synergy_search.sh`

and no command line arguments are specified. Successful completion of the script will be accompanied by the creation of numerous files, i.e.:

`ls data_bliss`

```
ALMANAC.CCLE_786O.csv ALMANAC.CCLE_A498.csv ALMANAC.CCLE_A549.csv ALMANAC.CCLE_ACHN.csv ...
```

and 

```
Merck.CCLE_A2058.csv Merck.CCLE_A2780.csv Merck.CCLE_A375.csv Merck.CCLE_A427.csv ...
```

where the prefix (i.e. `Merck` or `ALMANAC`) indicates the source data set and the suffix (i.e. `CCLE_786O`) indicates the particular cell line.

# Build the synergy prediction program

Both the single drug-based synergy prediction algorithm and the drug pair-based synergy prediction algorithm are implemented in a single C++ program, `predict_synergy`. This program requires the Message Passing Interface (MPI), and optionally OpenMP, for training and testing models on a computer cluster. The program has been developed and testing using the [OpenMPI](https://www.open-mpi.org/) implementation of MPI.

To build the `predict_synergy` program, run the `make` command in the `jdacs4c-staging/LANL_Synergy/predict` directory.

# Run the synergy prediction program

As mentioned above, the `predict_synergy` program uses MPI to train and test random forest-based predictive models of drug pair synergy. This requires running the program via the `mpirun` command that is included in all MPI distributions (however, some cluster scheduling systems may have other ways of launching MPI-based programs).

Running the `./predict_synergy` program with out any command line arguments lists all of the valid command line options:

```
Usage: predict_synergy version 0.6:
		[-o <output file>] (default is stdout)
		[--drug <file of drug features>]
		[--drug.random (randomize the drug features)]
		[--drug.blind (specific test set of drug features)]
		[--drug.target (per-drug target/mechanism of action)]
		[--cell <file of cell features>]
		[--cell.random (randomize the cell features)]
		[--cell.blind (specific test set of cell features)]
		[--cell.normalize (normalize the cell features to zero mean and unit variance)]
		[--enable-pair-prediction (perform both single drug and drug *pair* prediction)]
	Optional one hot encoding
		[--hot.cell (encode cell lines as one-hot features)]
		[--hot.tissue (encode cell line tissue type as one-hot features)]
		[--hot.drug (encode drugs as one-hot features)]
	Synergy variable to regress
		--synergy.train | --synergy <"cell name->file of training synergy measurements"> (can be repeated)
		[--synergy.test <"cell name->file of testing synergy measurements"> (can be repeated)]
		--synergy.average | --synergy.median | --synergy.fmean | --synergy.threshold <threshold value>
		[--synergy.NSC_to_CID (convert drug IDs from NSC to CID)]
		[--synergy.test.random (randomize the testing synergy variables to compute permutation p-value)]
		[--synergy.train.random (randomize the training synergy variables)]
		[-s <random number seed> (default is time based)]
	Cross-validation
		[--cv.folds <number of cross validation folds>] (default is 5)
	Generalization (when test and training are from different data sets)
		[--overlap.to_train] (shared cell lines/drugs assigend to train)
		[--overlap.to_test] (shared cell lines/drugs assigend to test)
	Regression parameters
		[--forest.size <number of trees>] (default is 100000)
		[--forest.bag <fraction of variables to bag>] (default is 0.3)
		[--forest.leaf <minimum leaf size>] (default is 3)
```

## Cross-validation on the NCI-ALMANAC data

An example of using `predict_synergy` to perform five-fold cross validation on the NCI-ALMANAC data set is provided in the 
shell script file `jdacs4c-staging/LANL_synergy/predict/grid_ALMANAC.sh`. Due to the length of time required to train the drug pair-based synergy models, this script was designed to run on a cluster computer made up of 18 servers, with 48, 1.9 GHz cores per server (for a total of 864 cores). Any cluster with a similar number of cores will take approximately *4 days* to generate the 25 random cross-validation samples for the following feature vector combinations:

* Full cell line and drug features
* Cell line features only
* Categorical tissue-type cell line features and drug features
* Categorical tissue-type cell line features only
* Drug features only
* No cell features and no drug features (as a negative control)

The `jdacs4c-staging/LANL_Synergy/predict/grid_ALMANAC.sh` file will need to be edited to adapt it to the local computing environment and to select between the Bliss and the Loewe synergy data.

The `predict_synergy` produces a rather verbose output file (with a summary of the input data and information about the results of each cross-validation fold). However, the most important information is presented at the very end of the output file. This information includes the Gini coefficients for both the drug pair-based (`Final drug pair-based synergy Gini`) and 
single drug-based (`Final single drug-based synergy Gini`) synergy models. The additional output information includes the "Area Under the Receiver Operator Curve" (AUROC) and the "enrichment" (a performance metric that is not used in the paper, since it is very similar to the Gini coefficient).

The Mean Squared Error (i.e. `MSE`) information in the output file reports on an effort to directly predict the single drug synergy probability (as a regression problem, rather than the classification problem of predicting drug pair synergy that is presented in the paper). However, this regression attempt was largely unsuccessful (although the predicted single drug probabilities do appear to be correlated with the actual single drug probabilities).

## Training on NCI-ALMANAC and testing on Merck

An example of using the `predict_synergy` program to train on the NCI-ALMANAC synergy data and test on the Merck synergy data is provided by the `jdacs4c-staging/LANL_Synergy/predict/batch_ALMANAC_to_Merck.sh` shell script file. Unlike the cross-validation script above, this script is much faster to execute (less than a day on the cluster computer detailed above), as it does not perform multiple iterations of cross validation.

By default, this script uses all of the specified cell lines for testing and training. However, as mentioned in the publication, there are some cell lines that are common to both the NCI-ALMANAC and the Merck data sets. To assign these shared cell lines to the NCI-ALMANAC training set, edit the `batch_ALMANAC_to_Merck.sh` script to add the `--overlap.to_train` flag. Alternatively, to assign these shared cell lines to the Merck testing set, edit the `batch_ALMANAC_to_Merck.sh` script to add the `--overlap.to_test` flag.


