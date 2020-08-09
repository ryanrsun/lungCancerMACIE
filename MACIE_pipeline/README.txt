The MACIE_pipeline folder provides the files necessary to annotate SNPs according to the trained model as described in our manuscript.

Included are:
- A snippet of the CADD v1.3 chromosome 22 file. To annotate arbitrary variants, you will need the full files for each chromosome, 
which can be downloaded from the CADD website. The full files are too big, so we have only provided the snippet needed for the reproducible
example.
- A class file noting the model specifications (which features to use, type of outcome, which classes they fall in). We provide this in the
/class folder.
- Various parameters files (provided in the /parameters folder) with the fitted model parameters and precalculated summary data measures.
You can also run the code in the "train" mode to retrain the model and estimate new model parameters.
- An example list of SNPs to be annotated, provided in the /snpList folder.
- An R file prepList_updated.R, which should be run (after you modify the appropriate input arguments at the top) to prepare the list of SNPs
for annotation. The R script will produce some intermediate files (annotating your SNPs, printing any SNPs with missing annotations) in the 
chosen output folder.
- The MACIE_iter200.py file, which was the code used to train the model and provide predictions.
- The file annotate.sh, which you can run after running prepList_updated.R (and modifying arguments) and will provide results.

In full, the procedure to reproduce our predictions is:
1. Create a list of SNPs (please use hg19) to be annotated and put it in the /snpLists folder (example provided).
2. Change file location inputs as described at top of prepList_updated.R (current arguments will work if you keep this same file structure) and run it.
3. Change file location inputs in annotate.sh and run it.
 
Note 1: The predictions header will look like "Chrom	Pos	Ref	Alt	(0, 1)	(1, 0)	(0, 0)	(1, 1)", where the first 0/1 refers to not-conserved/conserved and the second 0/1 refers to not-regulatory/regulatory (this can also be found based on the order of the class file). Therefore, (0, 1) refers to "not-conserved but regulatory"; (1, 0) refers to "conserved but not-regulatory"; (0, 0) refers to "not-conserved nor regulatory"; (1, 1) refers to "conserved and regulatory".
Note 2: The version of python used for this script is 2.7.14, later versions of python (2.7) should also work for this pipeline.
