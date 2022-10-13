
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

#

# Load libraries and dependencies ----------------------------------------------------

# R libraries by L.A. Courtenay

library(pValueRobust) # for p-Value calibrations
library(GraphGMM) # for Geometric Morphometrics

# to install these libraries use the install_github function from the devtools library
# library(devtools)
# install_github("LACourtenay/GraphGMM")
# install_github("LACourtenay/pValueRobust")

# external dependencies

library(geomorph) # for basic geometric morphometric functions
library(shapes) # for basic geometric morphometric functions and the o'higgins and dryden dataset
library(Morpho) # for the implementation of CVA and bgPCA
library(ggplot2) # for plotting
library(gridExtra) # for plotting
library(abind) # for matrix and tensor operations
library(caret) # for classification
library(gap) # for the CHOW test comparing regression lines
library(RVAideMemoire) # for MANOVA
library(umap) # for UMAP

#

# Load functions -----------------------------------------------------------

# function for the rotation of a point around a centroid.
# x is the x-coordinate of the point to be rotated
# y is the y-coordinate of the point to be rotated
# theta is the angle of rotation
# xa is the x-coordinate of the centroid
# ya is the y-coordinate of the centroid

rotate_point <- function(x, y, theta, xa = 0, ya = 0) {
  
  theta <- theta * (pi / 180)
  
  origin <- matrix(
    c(xa, ya),
    nrow = 2
  )
  
  rotation_matrix <- matrix(
    c(cos(theta), -sin(theta), sin(theta), cos(theta)),
    ncol = 2, byrow = TRUE
  )
  
  translated_to_origin <- matrix(
    c(x - xa, y - ya),
    nrow = 2
  )
  
  rotated_point <- origin + (rotation_matrix %*% translated_to_origin)
  
  return(t(rotated_point))
  
} 

# function for the creation of a base geometry
# n_landmarks are the number of landmarks required in the base geometry
# starting_point is the location of the first point that will be rotated to create the geometry
# NOTE: starting_point can be levaraged to control the centroid size of the geometry

create_base_shape <- function(n_landmarks, starting_point = c(0, 3)) {
  
  theta <- 360/n_landmarks
  
  base_shape <- array(dim = c(0, 2))
  base_shape <- abind::abind(base_shape, starting_point, along = 1)
  
  for (lm in 2:n_landmarks) {
    
    base_shape <- abind::abind(
      base_shape,
      rotate_point(
        base_shape[lm - 1, 1],
        base_shape[lm - 1, 2],
        theta
      ),
      along = 1
    )
    
  }
  
  return(base_shape)
  
}

# function for the shearing (deformation) of a shape
# shape is a matrix of coordinates that define the base shape
# epsilon is the degree of shearing
# axis is a integer value of either 1 or 2 that will shear the shape along the x (1) or y (2) axis

shear_shape <- function(shape, epsilon, axis) {
  
  if (axis == 1) {
    s_axis = 2
  } else {
    s_axis = 1
  }
  
  sheared_shape <- shape
  sheared_shape[,axis] <- shape[,axis] + (epsilon * shape[,s_axis])
  
  return(sheared_shape)
  
}

# function for the creation of a single sample of shapes based on the pertubation of landmarks for a single reference shape
# shape is a matrix of coordinate values defining the shape that will be used as reference to create the sample
# sigma is the deviation parameter of the gaussian distribution that is used for landmark deformation
# sample_size is the number of individuals to be created for this sample.
# NOTE: we recommend setting sigma to 1

create_single_sample <- function(shape, sigma, sample_size) {
  
  p <- nrow(shape)
  k <- ncol(shape)
  
  single_sample <- array(dim = c(p, k, 0))
  
  for (individual in 1:sample_size) {
    
    x_var <- c()
    y_var <- c()
    
    for (lm in 1:p) {
      
      x_var <- c(x_var, rnorm(1, mean = shape[lm,1], sd = sigma))
      y_var <- c(y_var, rnorm(1, mean = shape[lm,2], sd = sigma))
      
    }
    
    ind_var <- as.matrix(cbind(x_var, y_var), ncol = 2)
    
    single_individual <- shape + ind_var
    
    single_sample <- abind::abind(single_sample, single_individual, along = 3)
    
  }
  
  return(single_sample)
  
}

# a simple function for the calculation of a sample imbalance index based on a vector of factor labels
# this function calculates an adaptation of Shannon's entropy theory for the calculation of the balance of a dataset
# labels is a vector of factor labels

label_imbalance <- function(labels) {
  
  n <- length(labels)
  g <- length(levels(labels))
  c <- c()
  entropy_index <- c()
  for (ci in 1:g) {
    c <- c(
      c, table(labels)[ci][[1]]
    )
    entropy_index <- c(
      entropy_index, (c[ci] / n) * log(c[ci] / n)
    )
  }
  H <- -sum(entropy_index)
  balance_index <- H / log(g)
  return(balance_index)
  
}

# a simple adaptation of the previous function for the calculation of a sample imbalance index based on a vector of factor labels
# sample_sizes is a vector of sample sizes.

sample_imbalance <- function(sample_sizes) {
  
  n <- sum(sample_sizes)
  g <- length(sample_sizes)
  c <- sample_sizes
  entropy_index <- c()
  for (ci in 1:g) {
    entropy_index <- c(
      entropy_index, (c[ci] / n) * log(c[ci] / n)
    )
  }
  H <- -sum(entropy_index)
  balance_index <- H / log(g)
  return(balance_index)
  
}

# a function to create the experimental dataset with three sets of shapes for geometric morphometric analysis 
# shape1 is a matrix containing the base geometry that will be used to create the first sample
# shape2 is a matrix containing the base geometry that will be used to create the second sample
# shape3 is a matrix containing the base geometry that will be used to create the third sample
# epsilon defines the degree of shearing that will be performed on shape2 and shape3
# sigma is the deviation parameter for the gaussian pertubation of landmarks
# sample_sizes is a vector of length 3 that defines the size of each of the samples to be created

create_experimental_dataset <- function(shape1, shape2, shape3,
                                        epsilon, sigma = 1,
                                        sample_sizes = c(30, 30, 30)) {

  if (length(sample_sizes) != 3) {
    stop("the sample_size vector must be of length 3")
  }
  
  shape_prima_1 <- shear_shape(shape2, epsilon, 2)
  shape_prima_2 <- shear_shape(shape3, -epsilon, 2)
  
  group_1 <- create_single_sample(shape1, sigma, sample_sizes[1])
  group_2 <- create_single_sample(shape_prima_1, sigma, sample_sizes[2])
  group_3 <- create_single_sample(shape_prima_2, sigma, sample_sizes[3])
  
  experimental_sample <- abind::abind(group_1, group_2, group_3, along = 3)
  
  experimental_labels <- as.factor(c(
    rep("S1", dim(group_1)[3]),
    rep("S2", dim(group_2)[3]),
    rep("S3", dim(group_3)[3])
  ))
  
  return(list(
    coords = experimental_sample,
    labels = experimental_labels
  ))
  
}

# a function used to caclulate what the theoretical procrustes distances should be for a sample of sufficient sample size
# n_lm is the number of landmarks for the base geometry
# epsilon defines the degree of shearing that will be performed on samples 2 and 3
# beta is the parameter of Delta beta that controls the change in centroid size for analyses in form
# sigma is the deviation parameter for the gaussian pertubation of landmarks
# binary is an experimental parameter that can be used to condition experiments to only compare two reference shapes

theoretical_proc_d <- function(n_lm, epsilon, beta,
                               sigma = 1, binary = FALSE) {
  
  centroid_sizes = c(3, 3 + beta, 3 + (beta * 2))
  
  if (binary == TRUE) {
    
    shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
    shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
    
    dataset <- create_binary_dataset(shape1, shape2,
                                     epsilon = epsilon,
                                     sigma = sigma,
                                     sample_sizes = c(n_lm * 2,
                                                      n_lm * 2,
                                                      n_lm * 2))
    
  } else {
    
    shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
    shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
    shape3 <- create_base_shape(n_lm, c(0, centroid_sizes[3]))
    
    dataset <- create_experimental_dataset(shape1, shape2, shape3,
                                           epsilon = epsilon,
                                           sigma = sigma,
                                           sample_sizes = c(n_lm * 2,
                                                            n_lm * 2,
                                                            n_lm * 2))
    
  }
  
  cat("\nComputing theoretical procrustes distances...\n")
  
  if (beta > 0) {
    GPAform <- GraphGMM::GPA(dataset$coords, scale = FALSE)
  } else {
    GPAform <- GraphGMM::GPA(dataset$coords, scale = TRUE)
  }
  
  gdf <- geomorph::geomorph.data.frame(coords = GPAform$coordinates,
                                       Sample = dataset$labels)
  
  proc_dist_theo <- geomorph::procD.lm(coords ~ Sample, data = gdf)
  
  return(proc_dist_theo)
  
}

# the main function used to perform experiments for this study
# n_lm is the number of landmarks for the base geometry
# epsilon defines the degree of shearing that will be performed on samples 2 and 3
# beta is the parameter of Delta beta that controls the change in centroid size for analyses in form
# sample_sizes is a vector of length 3 that defines the size of each of the samples to be created
# sigma is the deviation parameter for the gaussian pertubation of landmarks
# cv is a boolean parameter that defines whether cross-validation is required for CVA and bgPCA
# pcs is a boolean parameter that defines whether PC scores are to be used as input for CVA and bgPCA
# centroid_plot is a boolean parameter that defines whether the entire scatter plot is to be visualised
#       or only the centroid of distributions

perform_experiment <- function(n_lm, epsilon, beta,
                               sample_sizes, sigma = 1,
                               cv = TRUE, pcs = TRUE,
                               centroid_plot = FALSE) {
  
  set.seed(NULL)
  # a function from Morpho (i think) seems to set the seed
  # therefore i reset the seed here.
  
  centroid_sizes = c(3, 3 + beta, 3 + (beta * 2))
  
  shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
  shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
  shape3 <- create_base_shape(n_lm, c(0, centroid_sizes[3]))
  
  cat("\nGenerating Samples")
  
  dataset <- create_experimental_dataset(shape1, shape2, shape3,
                                         epsilon = epsilon,
                                         sigma = sigma,
                                         sample_sizes = sample_sizes)
  
  sample_balance = round(label_imbalance(dataset$labels), 2)
  
  if (beta > 0) {
    cat("\nPerforming GPA in Form Space")
    GPAform <- GraphGMM::GPA(dataset$coords, scale = FALSE)
  } else {
    cat("\nPerforming GPA in Shape Space")
    GPAform <- GraphGMM::GPA(dataset$coords, scale = TRUE)
  }
  
  cat("\nCalculating PCA")
  
  pc_scores <- GraphGMM::pca_plot(
    GraphGMM::vector_from_landmarks(GPAform$coordinates)
  )$pc_scores
  
  expl_var <- GraphGMM::pca_plot(
    GraphGMM::vector_from_landmarks(GPAform$coordinates)
  )$variance
  
  eigenvals <- GraphGMM::pca_plot(
    GraphGMM::vector_from_landmarks(GPAform$coordinates)
  )$eigenvalues
  
  if (pcs == TRUE) {
    delta = round(
      ncol(pc_scores[,cumsum(expl_var) <= 0.95]) / min(sample_sizes), 2
    )
  } else {
    delta = round((n_lm*2) / min(sample_sizes), 2)
  }
  
  
  
  cat("\nCalculating CVA\n")
  
  if (pcs == TRUE) {
    if (cv == TRUE) {
      cva <- Morpho::CVA(pc_scores[,cumsum(expl_var) <= 0.95], dataset$labels,
                         cv = TRUE)
    } else {
      cva <- Morpho::CVA(pc_scores[,cumsum(expl_var) <= 0.95], dataset$labels,
                         cv = FALSE)
    }
    
  } else {
    if (cv == TRUE) {
      cva <- Morpho::CVA(GraphGMM::vector_from_landmarks(GPAform$coordinates),
                         dataset$labels,
                         cv = TRUE)
    } else {
      cva <- Morpho::CVA(GraphGMM::vector_from_landmarks(GPAform$coordinates),
                         dataset$labels,
                         cv = FALSE)
    }
    
  }
  
  cat("\nCalculating bgPCA")
  
  if (pcs == TRUE) {
    if (cv == TRUE) {
      bgpca <- Morpho::groupPCA(pc_scores[,cumsum(expl_var) <= 0.95],
                                dataset$labels,
                                cv = TRUE)
    } else {
      bgpca <- Morpho::groupPCA(pc_scores[,cumsum(expl_var) <= 0.95],
                                dataset$labels,
                                cv = FALSE)
    }
    
  } else {
    if (cv == TRUE) {
      bgpca <- Morpho::groupPCA(GraphGMM::vector_from_landmarks(GPAform$coordinates),
                                dataset$labels,
                                cv = TRUE)
    } else {
      bgpca <- Morpho::groupPCA(GraphGMM::vector_from_landmarks(GPAform$coordinates),
                                dataset$labels,
                                cv = FALSE)
    }
    
  }
  
  cat("\n\nEigenvalues and Canonical Roots: ")
  cat(paste0("\n\nPCA Eigenvalues: ", round(eigenvals[1], 2), " and ",
             round(eigenvals[2], 2)))
  cat(paste0("\nCVA Roots: ", round(cva$Var[1,1][[1]], 2), " and ",
             round(cva$Var[2,1][[1]], 2)))
  cat(paste0("\nbgPCA Eigenvalues: ", round(bgpca$Variance[1,1][[1]], 2), " and ",
             round(bgpca$Variance[2,1][[1]], 2)))
  
  cat("\n\nPlotting Results.")
  
  if (centroid_plot == TRUE) {
    
    gridExtra::grid.arrange(
      create_plot(GraphGMM::pca_plot(
        GraphGMM::vector_from_landmarks(GPAform$coordinates)
      ), dataset$labels, "pca", paste0(
        "\u03b5 = ", epsilon, ", \u0394\u03B2 = ", beta,
        ", max(\u03B4) = ", delta,
        ", B = ", sample_balance
      ), cv = FALSE, centroid_plot = TRUE),
      create_plot(cva, dataset$labels, "cva", "", cv = cv, centroid_plot = TRUE),
      create_plot(bgpca, dataset$labels, "bgpca", "", cv = cv, centroid_plot = TRUE),
      ncol = 3
    )
    
  } else {
    
    gridExtra::grid.arrange(
      GraphGMM::pca_plot(
        GraphGMM::vector_from_landmarks(GPAform$coordinates), dataset$labels,
        CI_ellipse = TRUE, point_size = 4, main = paste0(
          "\u03b5 = ", epsilon, ", \u0394\u03B2 = ", beta,
          ", max(\u03B4) = ", delta,
          ", B = ", sample_balance
        ))$pca_plot,
      create_plot(cva, dataset$labels, "cva", "", cv = cv, centroid_plot = FALSE),
      create_plot(bgpca, dataset$labels, "bgpca", "", cv = cv, centroid_plot = FALSE),
      ncol = 3
    )
    
  }
  
  proc_dist <- theoretical_proc_d(n_lm, epsilon, beta, sigma)
  
  cat("\n\nTheoretical Procrustes Distances:\n")
  
  cat(paste0("\nRsq: "), proc_dist$aov.table$Rsq[1])
  cat(paste0("\nF: "), proc_dist$aov.table$F[1])
  cat(paste0("\nZ: "), proc_dist$aov.table$Z[1])
  cat(paste0("\np: "), proc_dist$aov.table$`Pr(>F)`[1])
  
  cat("\n\nEmpirical Procrustes Distances:\n")
  
  gdf <- geomorph::geomorph.data.frame(coords = GPAform$coordinates,
                                       Sample = dataset$labels)
  
  proc_dist <- geomorph::procD.lm(coords ~ Sample, data = gdf)
  
  cat(paste0("\nRsq: "), proc_dist$aov.table$Rsq[1])
  cat(paste0("\nF: "), proc_dist$aov.table$F[1])
  cat(paste0("\nZ: "), proc_dist$aov.table$Z[1])
  cat(paste0("\np: "), proc_dist$aov.table$`Pr(>F)`[1])
  
  cat("\n\nComputing MANOVA on PCA...")
  
  pca_manova <- RVAideMemoire::pairwise.perm.manova(pc_scores[,1:2],
                                                    dataset$labels,
                                                    progress = FALSE,
                                                    test = "Wilks")
  
  cat("\nPCA MANOVA:\n")
  print(pca_manova$p.value)
  
  cat("\n\nComputing MANOVA on CVA...")
  
  if (cv == TRUE) {
    cva_manova <- RVAideMemoire::pairwise.perm.manova(cva$CVcv,
                                                      dataset$labels,
                                                      progress = FALSE,
                                                      test = "Wilks")
  } else {
    cva_manova <- RVAideMemoire::pairwise.perm.manova(cva$CVscores,
                                                      dataset$labels,
                                                      progress = FALSE,
                                                      test = "Wilks")
  }
  
  cat("\nCVA MANOVA:\n")
  print(cva_manova$p.value)
  
  cat("\n\nComputing MANOVA on bgPCA...")
  
  if (cv == TRUE) {
    bgpca_manova <- RVAideMemoire::pairwise.perm.manova(bgpca$CV,
                                                        dataset$labels,
                                                        progress = FALSE,
                                                        test = "Wilks")
  } else {
    bgpca_manova <- RVAideMemoire::pairwise.perm.manova(bgpca$Scores,
                                                        dataset$labels,
                                                        progress = FALSE,
                                                        test = "Wilks")
  }
  
  cat("\nbgPCA MANOVA:\n")
  print(bgpca_manova$p.value)
  
  cat("\n\nSimulation Done!\n")
  
}

# a function that can be used to create a plot of either a bgPCA, a CVA or a PCA
# data is the object produced by either the groupPCA, CVA or pca_plot functions.
# labels are the list of factor labels for the sample
# type is either "bgpca", "cva" or "pca"
# main is the title of the plot
# cv is a boolean value that indicates whether a Cv-bgPCA or Cv-CVA is to be visualised
# centroid_plot is a boolean parameter that defines whether the entire scatter plot is to be visualised
#       or only the centroid of distributions


create_plot <- function(data, labels, type, main, cv, centroid_plot) {
  
  if (type == "bgpca") {
    
    var1 <- data$Variance[1,2][[1]] * 100
    var2 <- data$Variance[2,2][[1]] * 100
    
    if (cv == TRUE) {
      data <- data$CV
      xlabel <- ggplot2::xlab(paste("Cv bgPC1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("Cv bgPC2 (",
                                    round(var2,2),"%)", sep = ""))
    } else {
      data <- data$Scores
      xlabel <- ggplot2::xlab(paste("bgPC1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("bgPC2 (",
                                    round(var2,2),"%)", sep = ""))
    }
    
    
  } else if (type == "cva") {
    
    var1 <- data$Var[1,2][[1]]
    var2 <- data$Var[2,2][[1]]
    
    if (cv == TRUE) {
      data <- data$CVcv
      xlabel <- ggplot2::xlab(paste("Cv CV1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("Cv CV2 (",
                                    round(var2,2),"%)", sep = ""))
    } else {
      data <- data$CVscores
      xlabel <- ggplot2::xlab(paste("CV1 (",
                                    round(var1,2),"%)", sep = ""))
      ylabel <- ggplot2::ylab(paste("CV2 (",
                                    round(var2,2),"%)", sep = ""))
    }
    
  } else {
    
    var1 <- data$variance[1] * 100
    var2 <- data$variance[2] * 100
    xlabel <- ggplot2::xlab(paste0("PC1 (", round(var1, 2), "%)"))
    ylabel <- ggplot2::ylab(paste0("PC2 (", round(var2, 2), "%)"))
    
    data <- data$pc_scores[,1:2]
    
  }
  
  data <- data.frame(
    x = data[,1],
    y = data[,2],
    Sample = as.factor(labels)
  )
  
  if (centroid_plot == FALSE) {
    base_plot <- ggplot2::ggplot(data = data,
                                 ggplot2::aes(x = x, y = y, colour = Sample))
    plot_colours <- ggplot2::scale_color_manual(values = c("black","red","blue",
                                                           "orange","darkgreen","violetred"))
    final_plot <- base_plot +
      ggplot2::geom_point(stat = "identity", size = 4) +
      plot_colours
    
    final_plot <- final_plot +
      ggplot2::stat_ellipse(size = 1) +
      xlabel + ylabel +
      ggplot2::ggtitle(main) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
        plot.title = ggplot2::element_text(face = "bold", size = 20),
        plot.subtitle = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA),
        axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
        axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
        legend.box.background = ggplot2::element_rect(colour = "black"),
        legend.position = "bottom"
      ) +
      ggplot2::geom_vline(xintercept = 0,
                          colour = "black",
                          size = 0.5,
                          linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0,
                          colour = "black",
                          linetype = "dashed",
                          size = 0.5)
  } else {
    
    centroids <- array(dim = c(3, 2)); target_index <- 1; for (target_sample in levels(data$Sample)) {
      centroids[target_index,] <- c(
        mean(data[data$Sample == target_sample, 1]), mean(data[data$Sample == target_sample, 2])
      )
      target_index <- target_index + 1
    }
    xlim <- c(min(data[,1]), max(data[,1]))
    ylim <- c(min(data[,1]), max(data[,1]))
    
    final_plot <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::geom_point(ggplot2::aes(x = centroids[1,1], y = centroids[1,2]), size = 7,
                          color = "black") +
      ggplot2::geom_point(ggplot2::aes(x = centroids[2,1], y = centroids[2,2]), size = 7,
                          color = "red") +
      ggplot2::geom_point(ggplot2::aes(x = centroids[3,1], y = centroids[3,2]), size = 7,
                          color = "blue") +
      xlabel + ylabel + ggplot2::ggtitle(main) +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
        plot.title = ggplot2::element_text(face = "bold", size = 20),
        plot.subtitle = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA),
        axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
        axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                             margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 18, face = "bold"),
        legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
        legend.box.background = ggplot2::element_rect(colour = "black"),
        legend.position = "bottom"
      ) +
      ggplot2::geom_vline(xintercept = 0,
                          colour = "black",
                          size = 0.5,
                          linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0,
                          colour = "black",
                          linetype = "dashed",
                          size = 0.5)
    
  }
  
  
  return(final_plot)
  
}

# a function that takes a 2x2 confusion matrix as input and calculates the precision

precision <- function(mat) {
  tp = mat[1, 1]
  fp = mat[1, 2]
  fn = mat[2, 1]
  tn = mat[2, 2]
  result = tp / (tp + fp)
  return(result)
}

# a function that takes a 2x2 confusion matrix as input and calculates the recall

recall <- function(mat) {
  tp = mat[1, 1]
  fp = mat[1, 2]
  fn = mat[2, 1]
  tn = mat[2, 2]
  result = tp / (tp + fn)
  return(result)
}

# a function that takes a 2x2 confusion matrix as input and calculates the f_stat

f_stat <- function(mat) {
  prec = precision(mat)
  rec = recall(mat)
  f = (2 * prec * rec) / (prec + rec)
  if(is.na(f)) {
    return(0)
  } else {
    return(f)
  }
}

# a function that takes a 2x2 confusion matrix as input and calculates the accuracy

accuracy <- function(mat) {
  tp = mat[1, 1]
  fp = mat[1, 2]
  fn = mat[2, 1]
  tn = mat[2, 2]
  result = (tp + tn) / (tp + fp + fn + tn)
  return(result)
}

# a function used to calculate correlation coeficients comparing two variables
# x is the x variable
# y is the y variable
# exponential is a boolean value that defines whether an exponential correlation model should be used.

correlation_calc <- function(x, y, exponential = FALSE) {
  
  if (exponential == FALSE) {
    
    if(shapiro.test(y)$p.value < 0.003) {
      
      cor_results <- suppressWarnings(cor.test(x, y, method = "kendall"))
      tau <- cor_results$estimate
      z <- cor_results$statistic
      p <- cor_results$p.value
      
      return(list(method = "Kendall", tau = tau, z = z, p = p))
      
    } else {
      
      cor_results <- suppressWarnings(cor.test(x, y, method = "pearson"))
      rho <- cor_results$estimate
      t <- cor_results$statistic
      p <- cor_results$p.value
      
      return(list(method = "Pearson", rho = rho, t = t, p = p))
      
    }
  } else {
    
    y[y == 0] <- 1e-05
    
    if(shapiro.test(y)$p.value < 0.003) {
      
      cor_results <- suppressWarnings(cor.test(x, log(y), method = "kendall"))
      tau <- cor_results$estimate
      z <- cor_results$statistic
      p <- cor_results$p.value
      
      return(list(method = "Kendall", tau = tau, z = z, p = p))
      
    } else {
      
      cor_results <- suppressWarnings(cor.test(x, log(y), method = "pearson"))
      rho <- cor_results$estimate
      t <- cor_results$statistic
      p <- cor_results$p.value
      
      return(list(method = "Pearson", rho = rho, t = t, p = p))
      
    }
  }
  
  
  
}

# a function used to calculate and print out the correlation between classification algorithm performance and balance of delta

calculate_class_eval_stats <- function(res_table) {
  
  cat("\nLDA vs B Correlation Results (Linear Model):\n")
  
  res <- correlation_calc(res_table$B, res_table$lda_f)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nPLSDA vs B Correlation Results (Linear Model):\n")
  
  res <- correlation_calc(res_table$B, res_table$plsda_f)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nLDA vs delta Correlation Results (Linear Model):\n")
  
  res <- correlation_calc(res_table$delta, res_table$lda_f)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho))
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau))
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nPLSDA vs delta Correlation Results (Linear Model):\n")
  
  res <- correlation_calc(res_table$delta, res_table$plsda_f)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho))
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau))
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nLDA vs B Correlation Results (Exponential Model):\n")
  
  res <- correlation_calc(res_table$B, res_table$lda_f, exponential = TRUE)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nPLSDA vs B Correlation Results (Exponential Model):\n")
  
  res <- correlation_calc(res_table$B, res_table$plsda_f, exponential = TRUE)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau * -1)) # times by -1 because the scales are reversed in delta
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nLDA vs delta Correlation Results (Exponential Model):\n")
  
  res <- correlation_calc(res_table$delta, res_table$lda_f, exponential = TRUE)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho))
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau))
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nPLSDA vs delta Correlation Results (Exponential Model):\n")
  
  res <- correlation_calc(res_table$delta, res_table$plsda_f, exponential = TRUE)
  
  if (res$method == "Pearson") {
    cat(paste0("\nrho: ", res$rho))
    cat(paste0("\nt: ", res$t))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  } else {
    cat(paste0("\nTau: ", res$tau))
    cat(paste0("\nz: ", res$z))
    cat(paste0("\np: ", res$p))
    cat(paste0("\nFPR: ", pValueRobust::FPR(res$p)))
  }
  
  cat("\n\nStatistical differences between LDA and PLSDA according to B:\n")
  
  chow_res <- gap::chow.test(
    res_table[,2], res_table[,4], res_table[,2], res_table[,6]
  )
  
  cat(paste0("\nChow's F: ", chow_res[1][[1]]))
  cat(paste0("\nChow's p: ", chow_res[4][[1]]))
  
  cat("\n\nStatistical differences between LDA and PLSDA according to delta:\n")
  
  chow_res <- gap::chow.test(
    res_table[,1], res_table[,4], res_table[,1], res_table[,6]
  )
  
  cat(paste0("\nChow's F: ", chow_res[1][[1]]))
  cat(paste0("\nChow's p: ", chow_res[4][[1]]))
  
  cat(paste0("\n\nMaximum LDA Accuracy:", max(res_table$lda_acc),"\n"))
  cat(paste0("Maximum LDA F-stat:", max(res_table$lda_f),"\n"))
  cat(paste0("\nMaximum PLSDA Accuracy:", max(res_table$plsda_acc),"\n"))
  cat(paste0("Maximum PLSDA F-stat:", max(res_table$lda_f),"\n"))
  
}

# a function used to take a 3 x 3 confusion matrix and convert it into a binary 2 x 2 confusion matrix
# for the study of evaluation metrics for the minority class of an imbalanced dataset

convert_to_binary_conf_mat <- function(conf_mat) {
  
  tp <- conf_mat[1,1]
  fp <- sum(conf_mat[1,2:3])
  fn <- sum(conf_mat[2:3,1])
  tn <- sum(conf_mat[2:3,2:3])
  
  return(
    matrix(
      c(tp, fp, fn, tn), byrow = TRUE, ncol = 2, nrow = 2
    )
  )
  
}

# the second most important function of this study, used to perform classification experiments
# n_lm is the number of landmarks for the base geometry
# epsilon defines the degree of shearing that will be performed on samples 2 and 3
# beta is the parameter of Delta beta that controls the change in centroid size for analyses in form
# sigma is the deviation parameter for the gaussian pertubation of landmarks
# delta_correct is a boolean value that defines whether delta values should always be small enough for a suitable application of LDA or PLSDA
# pcs is a boolean value that defines whether pc scores should be used as input to the classification algorithms

perform_classification_experiment <- function(n_lm, epsilon, beta, sigma = 1,
                                              delta_correct = FALSE,
                                              pcs = TRUE) {
  set.seed(NULL)
  
  centroid_sizes = c(3, 3 + beta, 3 + (beta * 2))
  
  sample_size_base <- (n_lm * 2) + (n_lm * 3)
  
  if (delta_correct == TRUE) {
    sample_size_base <- sample_size_base * 2
  }
  
  shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
  shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
  shape3 <- create_base_shape(n_lm, c(0, centroid_sizes[3]))
  
  trainset <- create_experimental_dataset(shape1, shape2, shape3,
                                          epsilon, sigma = sigma,
                                          sample_sizes = c(
                                            sample_size_base,
                                            sample_size_base,
                                            sample_size_base
                                          ))
  
  testset <- create_experimental_dataset(shape1, shape2, shape3,
                                         epsilon, sigma = sigma,
                                         sample_sizes = c(
                                           30,
                                           30,
                                           30
                                         ))
  
  all_data_set <- abind::abind(trainset$coords, testset$coords,
                               along = 3)
  
  if (delta_correct == TRUE) {
    sample_size_intervals <- rev(round(seq(sample_size_base / 2,
                                           sample_size_base,
                                           length.out = 20)))
  } else {
    sample_size_intervals <- rev(round(seq(ceiling(n_lm / 2),
                                           sample_size_base,
                                           length.out = 20)))
  }
  
  if (beta > 0) {
    cat("\nPerforming GPA in Form Space\n\n")
    GPAdata <- GraphGMM::GPA(all_data_set, scale = FALSE)
  } else {
    cat("\nPerforming GPA in Shape Space\n\n")
    GPAdata <- GraphGMM::GPA(all_data_set)
  }
  
  testset$coords <- GPAdata$coordinates[
    ,,(dim(all_data_set)[3] - 89):dim(all_data_set)[3]
  ]
  trainset$coords <- GPAdata$coordinates[
    ,,1:(dim(all_data_set)[3] - 90)
  ]
  
  train_vector <- GraphGMM::vector_from_landmarks(trainset$coords)
  test_vector <- GraphGMM::vector_from_landmarks(testset$coords)
  
  if (pcs == TRUE) {
    
    pca_results <- prcomp(as.matrix(train_vector), center = TRUE, scale = FALSE)
    
    expl_var <- pca_results$sdev^2 / sum(pca_results$sdev^2)
    test_pca <- data.frame(predict(pca_results, test_vector)[,cumsum(expl_var) <= 0.95],
                           Sample = as.factor(testset$labels))
    train_pca <- data.frame(pca_results$x[,cumsum(expl_var) <= 0.95],
                            Sample = as.factor(trainset$labels))
    
    n_var <- ncol(pca_results$x[,cumsum(expl_var) <= 0.95])
    
  } else {
    
    test_pca <- data.frame(test_vector,
                           Sample = as.factor(testset$labels))
    train_pca <- data.frame(train_vector,
                            Sample = as.factor(trainset$labels))
    
    n_var <- n_lm * 2
    
  }
  
  delta = c()
  B = c()
  lda_acc = c()
  lda_f = c()
  plsda_acc = c()
  plsda_f = c()
  
  cat("\nPerforming experiment:\n\n")
  
  pb <- txtProgressBar(min = 0, max = length(sample_size_intervals),
                       style = 3, width = 100, char = "=")
  prog = 1
  
  for (experiment in sample_size_intervals) {
    
    s1 <- train_pca[train_pca$Sample == "S1",]
    s2 <- train_pca[train_pca$Sample != "S1",]
    
    s1 <- s1[sample(nrow(s1), experiment, replace = FALSE),]
    
    dataset <- rbind(s1, s2)
    
    delta = c(delta, round((n_var) / nrow(s1), 2))
    B = c(B, label_imbalance(dataset$Sample))
    
    lda_model = suppressWarnings(
      caret::train(Sample ~ ., data = dataset, method="lda",
                   trControl = caret::trainControl(method = "cv"))
    )
    
    lda_results <- convert_to_binary_conf_mat(table(
      predict(lda_model, test_pca), test_pca$Sample
    ))
    
    lda_acc <- c(lda_acc, accuracy(lda_results))
    lda_f <- c(lda_f, f_stat(lda_results))
    
    plsda_model = plsda(dataset[,1:ncol(dataset)-1], dataset$Sample,
                        trControl = caret::trainControl(method = "cv"))
    
    plsda_results <- convert_to_binary_conf_mat(table(
      predict(plsda_model, test_pca[,1:ncol(test_pca)-1]), test_pca$Sample
    ))
    
    plsda_acc <- c(plsda_acc, accuracy(plsda_results))
    plsda_f <- c(plsda_f, f_stat(plsda_results))
    
    setTxtProgressBar(pb, prog)
    prog = prog + 1
    
  }
  
  result_table <- data.frame(
    delta = delta,
    B = B,
    lda_acc = lda_acc,
    lda_f = lda_f,
    plsda_acc = plsda_acc,
    plsda_f = plsda_f
  )
  
  colour_legend <- c("LDA" = "black", "PLSDA" = "red")
  bf_plot <- ggplot2::ggplot(data = result_table) +
    ggplot2::geom_point(ggplot2::aes(x = B, y = lda_f, color = "LDA"),
                        size = 4) +
    ggplot2::geom_point(ggplot2::aes(x = B, y = plsda_f,
                                     color = "PLSDA"),
                        size = 4) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = B, y = lda_f,
                                      color = "LDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = B, y = plsda_f,
                                      color = "PLSDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::scale_color_manual(values = colour_legend) +
    ggplot2::xlab("B") +
    ggplot2::ylab("F Measure") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 10, r = 0, b = 5, l = 0
                                           )),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 0, r = 10, b = 0, l = 0
                                           )),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    ) +
    ggplot2::scale_x_reverse()
  
  df_plot <- ggplot2::ggplot(data = result_table) +
    ggplot2::geom_point(ggplot2::aes(x = delta, y = lda_f, color = "LDA"),
                        size = 4) +
    ggplot2::geom_point(ggplot2::aes(x = delta, y = plsda_f,
                                     color = "PLSDA"),
                        size = 4) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = delta, y = lda_f,
                                      color = "LDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = delta, y = plsda_f,
                                      color = "PLSDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::scale_color_manual(values = colour_legend) +
    ggplot2::xlab("\u03B4") +
    ggplot2::ylab("F Measure") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 10, r = 0, b = 5, l = 0
                                           )),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 0, r = 10, b = 0, l = 0
                                           )),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    )
  
  calculate_class_eval_stats(result_table)
  
  cat("\n\nSimulation finished!\n\n")
  
  cat("\n\nPlotting results...\n")
  
  return(list(df_plot, bf_plot))
  
  #return(list(pca, bf_plot, result_table))
}

# this is a function for experimenting with tSNE and UMAP visualisations
# n_lm is the number of landmarks for the base geometry
# epsilon defines the degree of shearing that will be performed on samples 2 and 3
# beta is the parameter of Delta beta that controls the change in centroid size for analyses in form
# sigma is the deviation parameter for the gaussian pertubation of landmarks
# sample_sizes is a vector of length 3 that defines the size of each of the samples to be created
# pcs is a boolean value that defines whether pc scores should be used as input to the classification algorithms

tsne_umap <- function(n_lm, epsilon, beta,
                      sample_sizes, sigma = 1,
                      pcs = TRUE) {
  
  set.seed(NULL)
  # a function from Morpho (i think) seems to set the seed
  # therefore i reset the seed here.
  
  centroid_sizes = c(3, 3 + beta, 3 + (beta * 2))
  
  delta = round((n_lm*2) / min(sample_sizes), 2)
  
  shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
  shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
  shape3 <- create_base_shape(n_lm, c(0, centroid_sizes[3]))
  
  cat("\nGenerating Samples")
  
  dataset <- create_experimental_dataset(shape1, shape2, shape3,
                                         epsilon = epsilon,
                                         sigma = sigma,
                                         sample_sizes = sample_sizes)
  
  sample_balance = round(label_imbalance(dataset$labels), 2)
  
  if (beta > 0) {
    cat("\nPerforming GPA in Form Space")
    GPAform <- GraphGMM::GPA(dataset$coords, scale = FALSE)
  } else {
    cat("\nPerforming GPA in Shape Space")
    GPAform <- GraphGMM::GPA(dataset$coords, scale = TRUE)
  }
  
  cat("\nCalculating PCA")
  
  pc_scores <- GraphGMM::pca_plot(
    GraphGMM::vector_from_landmarks(GPAform$coordinates)
  )$pc_scores
  
  expl_var <- GraphGMM::pca_plot(
    GraphGMM::vector_from_landmarks(GPAform$coordinates)
  )$variance
  
  cat("\nCalculating tSNE\n")
  
  if (pcs == TRUE) {
    
    tsne <- GraphGMM::tsne_plot(as.data.frame(pc_scores[,cumsum(expl_var) <= 0.95]), dataset$labels,
                                CI_ellipse = TRUE, point_size = 4, main = "tSNE")
    
  } else {
    
    tsne <- GraphGMM::tsne_plot(GraphGMM::vector_from_landmarks(GPAform$coordinates),
                                dataset$labels,
                                CI_ellipse = TRUE, point_size = 4, main = "tSNE")
  }
  
  cat("\nCalculating UMAP")
  
  if (pcs == TRUE) {
    
    umap_data <- umap::umap(pc_scores[,cumsum(expl_var) <= 0.95])
    
  } else {
    
    umap_data <- umap::umap(GraphGMM::vector_from_landmarks(GPAform$coordinates))
    
  }
  
  umap_plot_data <- data.frame(
    x = umap_data$layout[,1],
    y = umap_data$layout[,2],
    Sample = as.factor(dataset$labels)
  )
  
  umap_plot <- ggplot2::ggplot(data = umap_plot_data,
                               ggplot2::aes(x = x, y = y, colour = Sample)) +
    ggplot2::geom_point(stat = "identity", size = 4) +
    ggplot2::scale_color_manual(values = c("black", "red", "blue",
                                           "orange", "grey")) +
    ggplot2::stat_ellipse(size = 1) +
    ggplot2::xlab("x") + ggplot2::ylab("y") +
    ggplot2::ggtitle("UMAP") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 10, r = 0, b = 5, l = 0)),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    ) +
    ggplot2::geom_vline(xintercept = 0,
                        colour = "black",
                        size = 0.5,
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "black",
                        linetype = "dashed",
                        size = 0.5)
  
  cat("\nPlotting Results.")
  
  gridExtra::grid.arrange(
    GraphGMM::pca_plot(
      GraphGMM::vector_from_landmarks(GPAform$coordinates), dataset$labels,
      CI_ellipse = TRUE, point_size = 4, main = paste0(
        "\u03b5 = ", epsilon, ", \u0394\u03B2 = ", beta,
        ", max(\u03B4) = ", delta,
        ", B = ", sample_balance
      ))$pca_plot,
    tsne$tSNE_plot,
    umap_plot, ncol = 3
  )
  
  cat("\n\nSimulation Done!\n")
  
}

# this is an experimental function that was not used in the main study
# this function is more oriented towards comparing datasets with just two samples (g = 2)
# distribution is the vector of data to be visualised
# labels are the list of factor labels for the sample
# main is the title of the plot

density_plot <- function(distribution, labels, main) {
  
  sample_levels <- levels(labels)
  dist1 <- as.numeric(distribution[labels == sample_levels[1],1])
  dist2 <- as.numeric(distribution[labels == sample_levels[2],1])
  
  dist1_dense <- density(dist1)
  dist2_dense <- density(dist2)
  
  ylim = c(
    0, max(c(dist1_dense$y, dist2_dense$y))
  )
  xlim = c(
    min(c(dist1_dense$x, dist2_dense$x)),
    max(c(dist1_dense$x, dist2_dense$x))
  )
  
  
  
  plot(dist1_dense$x, dist1_dense$y,
       xlim = xlim, ylim = ylim,
       type = "l", lwd = 2,
       xlab = "", ylab = "",
       xaxt = "none", yaxt = "none",
       main = main, cex.main = 2)
  lines(dist2_dense$x, dist2_dense$y,
        lwd = 2, col =" red", lty = 1)
  
  mtext(side = 1, line = 3, "Distribution", cex = 1.25, font = 2)
  mtext(side = 2, line = 3, "Density", cex = 1.25, font = 2)
  axis(1, round(seq(xlim[1], xlim[2], length.out = 5), 2), font = 1, cex.axis = 1)
  axis(2, round(seq(ylim[1], ylim[2], length.out = 5), 2), font = 1, cex.axis = 1)
  
  legend("topright", legend = sample_levels,
         inset = 0.015,
         col = c("Black", "Red"),
         lty = c(1,1), lwd = 2)
  
}

# this is an adaptation of create_experimental_dataset but for the comparison of only two samples (g = 2)

create_binary_dataset <- function(shape1, shape2,
                                  epsilon, sigma,
                                  sample_sizes = c(30, 30)) {
  
  shape_prima <- shear_shape(shape2, epsilon, 2)
  
  group_1 <- create_single_sample(shape1, sigma, sample_sizes[1])
  group_2 <- create_single_sample(shape_prima, 1, sample_sizes[2])
  
  experimental_sample <- abind::abind(group_1, group_2, along = 3)
  
  experimental_labels <- as.factor(c(
    rep("S", dim(group_1)[3]),
    rep("S'", dim(group_2)[3])
  ))
  
  return(list(
    coords = experimental_sample,
    labels = experimental_labels
  ))
  
}

# this is an adaptation of perform_experiment but oriented towards the study of two samples (g = 2)

perform_binary_experiment <- function(n_lm, beta, epsilon, sample_sizes, sigma = 1) {
  
  set.seed(NULL)
  
  centroid_sizes = c(3, 3 + beta)
  
  delta = round((n_lm*2) / min(sample_sizes), 2)
  
  cat("\nGenerating Samples")
  
  shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
  shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
  
  binary_dataset <- create_binary_dataset(shape1, shape2, epsilon, sigma = sigma,
                                          sample_sizes = sample_sizes)
  
  sample_balance = round(label_imbalance(binary_dataset$labels), 2)
  
  if (beta > 0) {
    cat("\nPerforming GPA in Form Space")
    GPAdata <- GraphGMM::GPA(binary_dataset$coords, scale = FALSE)
  } else {
    cat("\nPerforming GPA in Shape Space")
    GPAdata <- GraphGMM::GPA(binary_dataset$coords)
  }
  
  cat("\nCalculating PCA")
  
  pc_data <- GraphGMM::pca_plot(GraphGMM::vector_from_landmarks(GPAdata$coordinates),
                                binary_dataset$labels,
                                CI_ellipse = TRUE, point_size = 4, main = paste0(
                                  "\u03b5 = ", epsilon, ", \u0394\u03B2 = ", beta,
                                  ", \u03B4 \u003E ", delta,
                                  ", B = ", sample_balance
                                ))
  
  pc_scores <- pc_data$pc_scores
  expl_var <- pc_data$variance
  
  cat("\nCalculating CVA\n")
  
  if (pcs == TRUE) {
    cva <-Morpho::CVA(pc_scores[,cumsum(expl_var) <= 0.95],
                      binary_dataset$labels, plot = FALSE,
                      cv = TRUE)
  } else {
    cva <-Morpho::CVA(GraphGMM::vector_from_landmarks(GPAdata$coordinates),
                      binary_dataset$labels, plot = FALSE,
                      cv = TRUE)
  }
  
  cat("\nCalculating bgPCA")
  
  if (pcs == TRUE) {
    bgpca <- Morpho::groupPCA(pc_scores[,cumsum(expl_var) <= 0.95],
                              binary_dataset$labels)
  } else {
    bgpca <- Morpho::groupPCA(GraphGMM::vector_from_landmarks(GPAdata$coordinates),
                              binary_dataset$labels)
  }
  
  cat("\nVisualising results")
  
  par(mfrow = c(1,3), mar = c(5.1, 5, 4.1, 2.))  
  density_plot(pc_scores, binary_dataset$labels, "PCA")
  
  if (cv == TRUE) {
    density_plot(cva$CVcv, binary_dataset$labels, "Cv CVA")
    density_plot(bgpca$CV, binary_dataset$labels, "Cv bgPCA")
  } else {
    density_plot(cva$CVscores, binary_dataset$labels, "CVA")
    density_plot(bgpca$Scores, binary_dataset$labels, "bgPCA")
  }
  
  par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  
  cat(paste0(
    "\nSimulation summary: \n", paste0(
      "\u03b5 = ", epsilon, ", \u0394\u03B2 = ", beta,
      ", \u03B4 \u003E ", delta,
      ", B = ", sample_balance
    ), "\n"
  ))
  
  proc_dist <- theoretical_proc_d(n_lm, epsilon, beta, sigma,
                                  binary = TRUE)
  
  cat("\n\nTheoretical Procrustes Distances:\n")
  
  cat(paste0("\nRsq: "), proc_dist$aov.table$Rsq[1])
  cat(paste0("\nF: "), proc_dist$aov.table$F[1])
  cat(paste0("\nZ: "), proc_dist$aov.table$Z[1])
  cat(paste0("\np: "), proc_dist$aov.table$`Pr(>F)`[1])
  
  cat("\n\nEmpirical Procrustes Distances:\n")
  
  gdf <- geomorph::geomorph.data.frame(coords = GPAdata$coordinates,
                                       Sample = binary_dataset$labels)
  
  proc_dist <- geomorph::procD.lm(coords ~ Sample, data = gdf)
  
  cat(paste0("\nRsq: "), proc_dist$aov.table$Rsq[1])
  cat(paste0("\nF: "), proc_dist$aov.table$F[1])
  cat(paste0("\nZ: "), proc_dist$aov.table$Z[1])
  cat(paste0("\np: "), proc_dist$aov.table$`Pr(>F)`[1])
  
  cat("\n\nComputing Kruskal-Wallis on PCA...")
  
  pca_kw <- kruskal.test(pc_scores[,1],
                         binary_dataset$labels)
  
  cat("\nPCA Kruskal-Wallis rank sum test:\n")
  cat(paste0("\nChi2 = ", pca_kw$statistic))
  cat(paste0("\np = ", pca_kw$p.value))
  
  cat("\n\nComputing Kruskal-Wallis on CVA...")
  
  if (cv == TRUE) {
    cva_kw <- kruskal.test(cva$CVcv,
                           binary_dataset$labels)
  } else {
    cva_kw <- kruskal.test(cva$CVscores,
                           binary_dataset$labels)
  }
  
  cat("\nCVA Kruskal-Wallis rank sum test:\n")
  cat(paste0("\nChi2 = ", cva_kw$statistic))
  cat(paste0("\np = ", cva_kw$p.value))
  
  cat("\n\nComputing Kruskal-Wallis on bgPCA...")
  
  if (cv == TRUE) {
    bgpca_kw <- kruskal.test(bgpca$CV,
                             binary_dataset$labels)
  } else {
    bgpca_kw <- kruskal.test(bgpca$Scores,
                             binary_dataset$labels)
  }
  
  cat("\nbgPCA Kruskal-Wallis rank sum test:\n")
  cat(paste0("\nChi2 = ", bgpca_kw$statistic))
  cat(paste0("\np = ", bgpca_kw$p.value))
  
  cat("\nSimulation Done!")
  
}

# this is an adaptation of perform_classification_experiment but oriented towards the study of two samples (g = 2)

perform_binary_classification_experiment <- function(n_lm, epsilon, beta, sigma = 1,
                                                     delta_correct = FALSE,
                                                     pcs = TRUE) {
  
  set.seed(NULL)
  
  centroid_sizes = c(3, 3 + beta)
  
  sample_size_base <- (n_lm * 2) + (n_lm * 3)
  
  if (delta_correct == TRUE) {
    sample_size_base <- sample_size_base * 2
  }
  
  shape1 <- create_base_shape(n_lm, c(0, centroid_sizes[1]))
  shape2 <- create_base_shape(n_lm, c(0, centroid_sizes[2]))
  
  trainset <- create_binary_dataset(shape1, shape2,
                                    epsilon, sigma = sigma,
                                    sample_sizes = c(
                                      sample_size_base,
                                      sample_size_base
                                    ))
  
  testset <- create_binary_dataset(shape1, shape2,
                                   epsilon, sigma = sigma,
                                   sample_sizes = c(
                                     30,
                                     30
                                   ))
  
  all_data_set <- abind::abind(trainset$coords, testset$coords,
                               along = 3)
  
  if (delta_correct == TRUE) {
    sample_size_intervals <- rev(round(seq(sample_size_base / 2,
                                           sample_size_base,
                                           length.out = 20)))
  } else {
    sample_size_intervals <- rev(round(seq(ceiling(n_lm / 2),
                                           sample_size_base,
                                           length.out = 20)))
  }
  
  
  
  if (beta > 0) {
    cat("\nPerforming GPA in Form Space\n\n")
    GPAdata <- GraphGMM::GPA(all_data_set, scale = FALSE)
  } else {
    cat("\nPerforming GPA in Shape Space\n\n")
    GPAdata <- GraphGMM::GPA(all_data_set)
  }
  
  testset$coords <- GPAdata$coordinates[
    ,,(dim(all_data_set)[3] - 59):dim(all_data_set)[3]
  ]
  trainset$coords <- GPAdata$coordinates[
    ,,1:(dim(all_data_set)[3] - 60)
  ]
  
  train_vector <- GraphGMM::vector_from_landmarks(trainset$coords)
  test_vector <- GraphGMM::vector_from_landmarks(testset$coords)
  
  if (pcs == TRUE) {
    
    pca_results <- prcomp(as.matrix(train_vector), center = TRUE, scale = FALSE)
    
    expl_var <- pca_results$sdev^2 / sum(pca_results$sdev^2)
    test_pca <- data.frame(predict(pca_results, test_vector)[,cumsum(expl_var) <= 0.95],
                           Sample = as.factor(testset$labels))
    train_pca <- data.frame(pca_results$x[,cumsum(expl_var) <= 0.95],
                            Sample = as.factor(trainset$labels))
  } else {
    
    test_pca <- data.frame(test_vector,
                           Sample = as.factor(testset$labels))
    train_pca <- data.frame(train_vector,
                            Sample = as.factor(trainset$labels))
    
  }
  
  delta = c()
  B = c()
  lda_acc = c()
  lda_f = c()
  plsda_acc = c()
  plsda_f = c()
  
  cat("\nPerforming experiment:\n\n")
  
  pb <- txtProgressBar(min = 0, max = length(sample_size_intervals),
                       style = 3, width = 100, char = "=")
  prog = 1
  
  for (experiment in sample_size_intervals) {
    
    s1 <- train_pca[train_pca$Sample == "S",]
    s2 <- train_pca[train_pca$Sample == "S'",]
    
    s1 <- s1[sample(nrow(s1), experiment, replace = FALSE),]
    
    dataset <- rbind(s1, s2)
    
    delta = c(delta, round((n_lm*2) / nrow(s1), 2))
    B = c(B, label_imbalance(dataset$Sample))
    
    lda_model = suppressWarnings(
      caret::train(Sample ~ ., data = dataset, method="lda",
                   trControl = caret::trainControl(method = "cv"))
    )
    
    lda_results <- table(
      predict(lda_model, test_pca), test_pca$Sample
    )
    
    lda_acc <- c(lda_acc, accuracy(lda_results))
    lda_f <- c(lda_f, f_stat(lda_results))
    
    plsda_model = plsda(dataset[,1:ncol(dataset)-1], dataset$Sample,
                        trControl = caret::trainControl(method = "cv"))
    
    plsda_results <- table(
      predict(plsda_model, test_pca[,1:ncol(test_pca)-1]), test_pca$Sample
    )
    
    plsda_acc <- c(plsda_acc, accuracy(plsda_results))
    plsda_f <- c(plsda_f, f_stat(plsda_results))
    
    setTxtProgressBar(pb, prog)
    prog = prog + 1
    
  }
  
  result_table <- data.frame(
    delta = delta,
    B = B,
    lda_acc = lda_acc,
    lda_f = lda_f,
    plsda_acc = plsda_acc,
    plsda_f = plsda_f
  )
  
  colour_legend <- c("LDA" = "black", "PLSDA" = "red")
  bf_plot <- ggplot2::ggplot(data = result_table) +
    ggplot2::geom_point(ggplot2::aes(x = B, y = lda_f, color = "LDA"),
                        size = 4) +
    ggplot2::geom_point(ggplot2::aes(x = B, y = plsda_f,
                                     color = "PLSDA"),
                        size = 4) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = B, y = lda_f,
                                      color = "LDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = B, y = plsda_f,
                                      color = "PLSDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::scale_color_manual(values = colour_legend) +
    ggplot2::xlab("B") +
    ggplot2::ylab("F Measure") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 10, r = 0, b = 5, l = 0
                                           )),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 0, r = 10, b = 0, l = 0
                                           )),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    ) +
    ggplot2::scale_x_reverse()
  
  df_plot <- ggplot2::ggplot(data = result_table) +
    ggplot2::geom_point(ggplot2::aes(x = delta, y = lda_f, color = "LDA"),
                        size = 4) +
    ggplot2::geom_point(ggplot2::aes(x = delta, y = plsda_f,
                                     color = "PLSDA"),
                        size = 4) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = delta, y = lda_f,
                                      color = "LDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::geom_smooth(method = "lm",
                         ggplot2::aes(x = delta, y = plsda_f,
                                      color = "PLSDA"),
                         se = FALSE,
                         formula = y ~ x) +
    ggplot2::scale_color_manual(values = colour_legend) +
    ggplot2::xlab("\u03B4") +
    ggplot2::ylab("F Measure") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 10, r = 0, b = 5, l = 0
                                           )),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18,
                                           margin = ggplot2::margin(
                                             t = 0, r = 10, b = 0, l = 0
                                           )),
      axis.text.x = ggplot2::element_text(angle = 90, size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 18, face = "bold"),
      legend.background = ggplot2::element_rect(fill = add_alpha("#CCCCCC", alpha = 0.2)),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.position = "bottom"
    )
  
  calculate_class_eval_stats(result_table)
  
  cat("\n\nSimulation finished!\n\n")
  
  cat("\n\nPlotting results...\n")
  
  return(list(df_plot, bf_plot))
  
  #return(list(pca, bf_plot, result_table))
  
}

# miscellaneous function for processing colours

add_alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

# a set of functions to calculate the precision, recall, accuracy and f1 scores for all samples in a confusion matrix.

evaluation_metrics <- function(confusion_matrix) {
  
  create_binary_table <- function(conf_mat, target) {
    
    tp <- conf_mat[target, target]
    fn <- sum(conf_mat[,target]) - tp
    fp <- sum(conf_mat[target,]) - tp
    tn <- sum(conf_mat) - sum(c(tp, fn, fp))
    
    return(
      matrix(
        c(tp, fp, fn, tn), byrow = TRUE, ncol = 2, nrow = 2
      )
    )
    
  }
  
  sample_names <- colnames(confusion_matrix)
  
  number_of_labels <- nrow(confusion_matrix)
  
  result_table <- data.frame(Sample = character(),
                             Accuracy = numeric(),
                             Precision = numeric(),
                             Recall = numeric(),
                             F1 = numeric())
  
  for (target_class in 1:number_of_labels) {
    
    target_matrix <- create_binary_table(confusion_matrix, target_class)
    
    results <- c(
      sample_names[target_class],
      accuracy(target_matrix),
      precision(target_matrix),
      recall(target_matrix),
      f_stat(target_matrix)
    )
    
    
    result_table[nrow(result_table) + 1,] <- results
    
  }
  
  return(result_table)
  
}

#