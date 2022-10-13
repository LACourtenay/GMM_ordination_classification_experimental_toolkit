
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

# please ensure that the functions.R file has been sourced prior to running code in this file

#

# code used to produce figure 13 and 14 from the original paper

perform_experiment(
    n_lm = 100,
    epsilon = 0.0,
    beta = 0.0,
    sample_sizes = c(30, 200, 200),
    cv = FALSE,
    pcs = FALSE
)

tsne_umap(
    n_lm = 100,
    epsilon = 0.0,
    beta = 0.0,
    sample_sizes = c(30, 200, 200),
    pcs = FALSE
)

perform_experiment(
    n_lm = 100,
    epsilon = 0.2,
    beta = 0.0,
    sample_sizes = c(30, 200, 200),
    cv = FALSE,
    pcs = FALSE
)

tsne_umap(
    n_lm = 100,
    epsilon = 0.2,
    beta = 0.0,
    sample_sizes = c(30, 200, 200),
    pcs = FALSE
)

#