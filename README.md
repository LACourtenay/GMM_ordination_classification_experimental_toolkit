# GMM_ordination_classification_experimental_toolkit

A set of functions and files with code to perform theoretical experiments comparing the comparance of different ordination and classification tests for geometric morphometric data. This code was originally developed to analyse sample imbalance.

-----------------------------------------------------------------------------------------------------------------

<i>
Author

Lloyd A. Courtenay

Email

ladc1995@gmail.com

ORCID

https://orcid.org/0000-0002-4810-2001

Current Afiliations:

Universidad de Salamanca [USAL]

</i>

---------------------------------------------------------------------------------------------------

This code was designed and prepared for the study by:
<b> Courtenay, L.A. (2023)
Can we restore balance to Geometric Morphometrics? A theoretical evaluation of how sample imbalance conditions ordination and classification.
Evoltuionary Biology. 50:90-110. https://doi.org/10.1007/s11692-022-09590-0</b>

TIDOP research group website: http://tidop.usal.es/

---------------------------------------------------------------------------------------------------

Prior to using any of the code, the functions.R file should first be sourced (while the user must also first check that all libraries are installed). Once the functions.R file has been sourced all code from other files should work fine. The functions.R file additionally contains some details about how to use each of the functions, while the code used to generate figures from the aforementioned study has been included in each of the separate files;

* <b>basic_experiments.R</b> contains all code used to perform the main study
* <b>o_higgins_and_dryden_demonstration.R</b> contains the code to perform experiments with the O'higgins and Dryden (1993) apes dataset available from the <i>shapes</i> R library.
* <b>tsne_and_umap_experiments.R</b> contains all code used to demonstrate the power of TSNE and UMAP as an alternative for ordination
* <b>smote_and_adasyn_experiments.R</b> contains all code used to demonstrate the power of SMOTE and ADASYN for the balancing of datasets.
  
---------------------------------------------------------------------------------------------------

Comments, questions, doubts, suggestions and corrections can all be directed to L. A. Courtenay at the email provided above.
