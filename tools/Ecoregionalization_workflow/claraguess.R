##30/04/2025
##Jean Le Cras
### Clustering with Clara algorithm with an option to determine the optimal number of clusters based on SIH index

#load libraries
library(cluster)
library(dplyr)
library(tidyverse)

#load arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("This tool needs at least one argument")
}

#load data
enviro_path <- args[1]
preds_path <- args[2]
taxa_path <- args[3]
type <- args[4]
k <- as.integer(args[5])
metric <- args[6]
samples <- as.integer(args[7])
env.data <- read.table(enviro_path, sep = "\t", header = TRUE, dec = ".", na.strings = "-9999")

data_split = str_split(preds_path, ",")
preds.data = NULL

for (i in 1:length(data_split[[1]])) {
  df <- read.table(data_split[[1]][i], dec=".", sep="\t", header=T, na.strings="NA")
  preds.data <- rbind(preds.data, df)
  remove(df)
}

names(preds.data) <- c("lat", "long", "pred", "taxon")

##List of modelled taxa used for clustering
taxa_list <- tv <- read.table(taxa_path, dec=".", sep=" ", header=F, na.strings = "NA") 
names(taxa_list) <- c("a")

preds.data <- preds.data[which(preds.data$taxon %in% taxa_list$a),]

#format data
abondance_mat <- matrix(preds.data$pred, nrow = nrow(env.data), ncol = nrow(preds.data) / nrow(env.data))
abondance_mat <- data.frame(abondance_mat)

message("abondance_mat", abondance_mat)
message("colnames(abondance_mat)", colnames(abondance_mat))
message("unique(preds.data$taxon)", unique(preds.data$taxon))
names(abondance_mat) <- unique(preds.data$taxon)

#select the clara model with the optimal number of clusters
model <- NULL

if (type == "auto") {
  sih_scores <- c()
  models <- list()
  
  for (i in 2:k) {
    models[[i]] <- clara(abondance_mat, i, metric = metric, samples = samples, stand = TRUE)
    sih_scores[i] <- models[[i]]$silinfo$avg.width
  }
  png("sih_scores.png")
  plot(2:k, sih_scores[2:k], type = "b", xlab = "Number of clusters", ylab = "SIH index")
  dev.off()
  
  best_k <- which.max(sih_scores[2:k]) + 1
  model <- models[[best_k]]
  k <- best_k
} else {
  model <- clara(abondance_mat, k, metric = metric, samples = samples, stand = TRUE)
}

#saving results
png("silhouette_plot.png")
plot(silhouette(model), main = paste("Silhouette plot for k =", k))
dev.off()

full.data <- cbind(preds.data[1:nrow(abondance_mat), 1:2], model$clustering)
names(full.data) <- c("lat", "long", "cluster")
full.data <- cbind(full.data, abondance_mat, env.data[, 3:ncol(env.data)])

write.table(full.data[1:3], file = "data_cluster.tabular", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(full.data, file = "clustered_taxas_env.tabular", quote = FALSE, sep = "\t", row.names = FALSE)