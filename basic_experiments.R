
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

# please ensure that the functions.R file has been sourced prior to running code in this file

#

# ordination experiments ------------------------

# 8 landmarks

# no change in shape or form

perform_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = TRUE,
    cv = FALSE
)

perform_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = TRUE,
    cv = TRUE
)

perform_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = FALSE,
    cv = FALSE
)

# slight change in shape

perform_experiment(
    n_lm = 8,
    epsilon = 0.1,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = TRUE,
    cv = FALSE
)

perform_experiment(
    n_lm = 8,
    epsilon = 0.1,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = TRUE,
    cv = TRUE
)

perform_experiment(
    n_lm = 8,
    epsilon = 0.1,
    beta = 0.0,
    sample_sizes = c(16, 50, 50),
    pcs = FALSE,
    cv = FALSE
)

# slight change in form

perform_experiment(
  n_lm = 8,
  epsilon = 0.0,
  beta = 0.15,
  sample_sizes = c(16, 50, 50),
  pcs = TRUE,
  cv = FALSE
)

perform_experiment(
  n_lm = 8,
  epsilon = 0.0,
  beta = 0.15,
  sample_sizes = c(16, 50, 50),
  pcs = TRUE,
  cv = TRUE
)

perform_experiment(
  n_lm = 8,
  epsilon = 0.0,
  beta = 0.15,
  sample_sizes = c(16, 50, 50),
  pcs = FALSE,
  cv = FALSE
)

# 100 landmarks

# no change in shape or form

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = FALSE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = TRUE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = FALSE,
  cv = FALSE
)

# slight change in shape

perform_experiment(
  n_lm = 100,
  epsilon = 0.05,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = FALSE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.05,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = TRUE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.05,
  beta = 0.0,
  sample_sizes = c(200, 500, 500),
  pcs = FALSE,
  cv = FALSE
)

# slight change in form

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.05,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = FALSE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.05,
  sample_sizes = c(200, 500, 500),
  pcs = TRUE,
  cv = TRUE
)

perform_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.05,
  sample_sizes = c(200, 500, 500),
  pcs = FALSE,
  cv = FALSE
)

#

# classification experiments --------------------

# 8 landmarks

# no change in shape or form

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = TRUE
)

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = FALSE
)

# slight change in shape

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.01,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = TRUE
)

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.01,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = FALSE
)

# extreme change in shape

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.2,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = TRUE
)

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.2,
    beta = 0.0,
    sigma = 1,
    delta_correct = TRUE,
    pcs = FALSE
)

# slight change in form

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.01,
    sigma = 1,
    delta_correct = TRUE,
    pcs = TRUE
)

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.01,
    sigma = 1,
    delta_correct = TRUE,
    pcs = FALSE
)

# extreme change in form

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.2,
    sigma = 1,
    delta_correct = TRUE,
    pcs = TRUE
)

perform_classification_experiment(
    n_lm = 8,
    epsilon = 0.0,
    beta = 0.2,
    sigma = 1,
    delta_correct = TRUE,
    pcs = FALSE
)

# 100 landmarks

# no change in shape or form

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = TRUE
)

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = FALSE
)

# slight change in shape

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.01,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = TRUE
)

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.01,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = FALSE
)

# extreme change in shape

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.1,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = TRUE
)

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.1,
  beta = 0.0,
  sigma = 1,
  delta_correct = TRUE,
  pcs = FALSE
)

# slight change in form

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.01,
  sigma = 1,
  delta_correct = TRUE,
  pcs = TRUE
)

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.01,
  sigma = 1,
  delta_correct = TRUE,
  pcs = FALSE
)

# extreme change in form

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.1,
  sigma = 1,
  delta_correct = TRUE,
  pcs = TRUE
)

perform_classification_experiment(
  n_lm = 100,
  epsilon = 0.0,
  beta = 0.1,
  sigma = 1,
  delta_correct = TRUE,
  pcs = FALSE
)

#