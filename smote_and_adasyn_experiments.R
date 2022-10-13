
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

#

# load some additional dependencies ---------------------------

library(performanceEstimation) # for SMOTE
library(smotefamily) # for ADASYN

#

# prepare experimental data -------------------------

trial <- create_single_sample(create_base_shape(100), 1, 200)
trial2 <- create_single_sample(create_base_shape(100), 1, 500)
trial3 <- create_single_sample(create_base_shape(100), 1, 500)

n_pc_scores <- 83

data <- abind::abind(trial, trial2, trial3, along = 3)
labels <- c(rep("S1", 200), rep("S2", 500), rep("S3", 500))

GPAshape <- GPA(data)

#

# basic PCA visualisation ------------------------------

pca_data <- pca_plot(vector_from_landmarks(GPAshape$coordinates), as.factor(labels),
                     CI_ellipse = TRUE, point_size = 4, main = "PCA")
pca_data$pca_plot

#

# perform SMOTE ------------------------------

pc_scores <- pca_data$pc_scores[,1:n_pc_scores]
data_set <- data.frame(pc_scores, Sample = as.factor(labels))

smote_data <- smote(Sample ~ ., data_set, k = 5, perc.over = 5)

smote_plot_data <- data.frame(
  x = smote_data[smote_data$Sample == "S1",1],
  y = smote_data[smote_data$Sample == "S1",2],
  Sample = smote_data[smote_data$Sample == "S1",]$Sample
)

original_data <- data.frame(
  x = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",1],
  y = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",2],
  Sample = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",]$Sample
)

plot(smote_data[smote_data$Sample == "S1",1:2], col = "red", pch = 19,
     xlab = "PC1", ylab = "PC2", cex = 0.75)
points(pc_scores[data_set$Sample == "S1",1:2], pch = 19, cex = 1.5)

smote_plot_data <- data.frame(
  PC1 = smote_data[smote_data$Sample == "S1",1],
  PC2 = smote_data[smote_data$Sample == "S1",2],
  Sample = "S1'"
)

original_s1 <- as.data.frame(pc_scores[data_set$Sample == "S1",1:2])
colnames(original_s1) <- c("PC1", "PC2")
original_s1$Sample <- "S1"

ggplot2::ggplot() +
  ggplot2::geom_point(
    data = smote_plot_data, ggplot2::aes(x = PC1, y = PC2),
    colour = "red"
  ) +
  ggplot2::geom_point(
    data = original_s1, ggplot2::aes(x = PC1, y = PC2),
    colour = "black", size = 4
  ) +
  ggplot2::theme_bw() +
  ggplot2::xlab(
    paste0("PC1 (", round(pca_data$variance[1] * 100, 2), "%)")
  ) + ggplot2::ylab(
    paste0("PC2 (", round(pca_data$variance[2] * 100, 2), "%)")
  ) +
  ggplot2::ggtitle("Balanced S1 using SMOTE") +
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

smote_plot_data <- rbind(smote_plot_data, original_data)

ggplot2::ggplot(data = smote_plot_data,
                ggplot2::aes(x = x, y = y, colour = Sample)) +
  ggplot2::geom_point(stat = "identity", size = 4) +
  ggplot2::scale_color_manual(values = c("black","red","blue","orange",
                                         "darkgreen","violetred")) +
  ggplot2::stat_ellipse(size = 1) +
  ggplot2::xlab(
    paste0("PC1 (", round(pca_data$variance[1] * 100, 2), "%)")
  ) + ggplot2::ylab(
    paste0("PC2 (", round(pca_data$variance[2] * 100, 2), "%)")
  ) +
  ggplot2::ggtitle("SMOTE PCA") +
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

#

# perform ADASYN ------------------------------

target_class <- as.character(data_set$Sample)
target_class[target_class != "S1"] <- "p" # need to set the majority classes to "p"
target_class[target_class == "S1"] <- "n" # need to set the minority classes to "n"

trial_data <- data.frame(pc_scores, target = target_class)

adasyn_data <- ADAS(trial_data[,1:n_pc_scores], trial_data$target, K = 5)$syn_data[,1:n_pc_scores]

adasyn_plot_data <- data.frame(
  x = adasyn_data[,1],
  y = adasyn_data[,2],
  Sample = rep("S1", nrow(adasyn_data))
)

original_data <- data.frame(
  x = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",1],
  y = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",2],
  Sample = data_set[data_set$Sample == "S2" | data_set$Sample == "S3",]$Sample
)

adasyn_plot_data <- rbind(adasyn_plot_data, original_data)

ggplot2::ggplot(data = adasyn_plot_data,
                ggplot2::aes(x = x, y = y, colour = Sample)) +
  ggplot2::geom_point(stat = "identity", size = 4) +
  ggplot2::scale_color_manual(values = c("black","red","blue","orange",
                                         "darkgreen","violetred")) +
  ggplot2::stat_ellipse(size = 1) +
  ggplot2::xlab(
    paste0("PC1 (", round(pca_data$variance[1] * 100, 2), "%)")
  ) + ggplot2::ylab(
    paste0("PC2 (", round(pca_data$variance[2] * 100, 2), "%)")
  ) +
  ggplot2::ggtitle("ADASYN PCA") +
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

#