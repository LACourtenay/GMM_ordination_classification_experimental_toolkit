
#

# Code written by Lloyd A. Courtenay #
# Lloyd A. Courtenay - ladc1995@gmail.com (Universidad de Salamanca [USAL]) #

# please ensure that the functions.R file has been sourced prior to running code in this file

#

# load and prepare data ---------------------------------

data(apes)

# perform GPA

GPAshape <- GPA(apes$x)

# generate binary labels

labels <- as.character(apes$group)
labels[apes$group == "gorf" | apes$group == "gorm"] <- "Gorilla"
labels[apes$group == "pongof" | apes$group == "pongom"] <- "Orangutan"
labels[apes$group == "panf" | apes$group == "panm"] <- "Chimpanzee"
labels <- as.factor(labels)

#

# perform base CVA and bgPCA -------------------------

# CVA and bgPCA with separation between species and sexes

apes_cva <- Morpho::CVA(GPAshape$coordinates, apes$group)
apes_bgpca <- Morpho::groupPCA(GPAshape$coordinates, apes$group)

# CVA and bgPCA with separation between species but no separation for sex

apes_cva2 <- Morpho::CVA(GPAshape$coordinates, labels)
apes_bgpca2 <- Morpho::groupPCA(GPAshape$coordinates, labels)

#

# visualise results --------------------------------

# with division between the sexes

pca_plot(vector_from_landmarks(GPAshape$coordinates), apes$group, CI_ellipse = TRUE,
         point_size = 4, main = "PCA")

create_plot(apes_cva, apes$group, "cva", "CVA", cv = FALSE)

create_plot(apes_bgpca, apes$group, "bgpca", "bgPCA", cv = FALSE)

# without the division between the sexes

pca_plot(vector_from_landmarks(GPAshape$coordinates), labels, CI_ellipse = TRUE,
         point_size = 4, main = "PCA")

create_plot(apes_cva2, labels, "cva", "CVA", cv = FALSE)
create_plot(apes_cva2, apes$group, "cva", "CVA", cv = FALSE)

create_plot(apes_bgpca2, labels, "bgpca", "bgPCA", cv = FALSE)
create_plot(apes_bgpca2, apes$group, "bgpca", "bgPCA", cv = FALSE)

#