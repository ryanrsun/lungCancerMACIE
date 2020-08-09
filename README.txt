This folder contains data used in "Integration of multi-omic annotation data to 
prioritize and characterize inflammation and immune-related risk variants in
 squamous cell lung cancer," Sun et al. (2020+).

Included are the MACIE predictions for all significant coding and noncoding SNPs
(see lung_cancer_coding_MACIE.txt and lung_cancer_noncoding_MACIE.txt). In these
files, the (0,1) column is the probability of the regulatory class only, the 
the (1,0) column is the probability of the conserved conserved class only, the 
(0,0) column is the probability of neither class, and the (1,1) column is the probability
of both classes. For coding SNPs, the regulatory class is replaced by the protein 
function class.

The full annotation data for all these SNPs can be found in lung_cancer_sig_annotations.txt.

The gene definitions we used can be found in hg19_gene_definitions.txt.

The code to reproduce our results or annotate your own SNPs is in the MACIE_pipeline folder.

Full MACIE predictions for all noncoding SNPs in the CADD database (approximately 7.8
million SNPs) are publicly available for lookup/download from dropbox at the links in
MACIE_dropbox_links.txt.

