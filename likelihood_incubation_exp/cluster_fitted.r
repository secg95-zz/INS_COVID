library(RJSONIO)
library(factoextra)

root = "likelihood_incubation_exp/tuned/"
model_subdirs = list.dirs(root, recursive=FALSE)
model_subdirs = model_subdirs[!grepl("legacy", model_subdirs)]

models = list()
# unpack each fitted model's parameters
for (model_subdir in model_subdirs) {
  model_path = paste(model_subdir, "model.json", sep="/")
  model = fromJSON(model_path)
  as_row = c(model$beta, model$tau1, model$tau2, model$N0, model$A0, model$loss)
  models[[model_subdir]] = as_row
}
models = do.call("rbind", models)
models = data.frame(models)
colnames(models) = c(paste0("beta", 1:length(model$beta)), "tau1", "tau2", "N0", "A0", "loss")
X = models[colnames(models)[colnames(models) != "loss"]]
col_sds = apply(X, 2, sd)
col_means = apply(X, 2, mean)
X = sweep(X, 2, col_means, `-`)
X = sweep(X, 2, col_sds, `/`)

# perform k-means clustering at several numbers of clusters
wss = c()
clusterings = list()
for (k in 1:20) {
  clusters = kmeans(X, k, nstart=1000)
  wss = c(wss, clusters$tot.withinss)
  # de-normalize clusters
  clusters$centers = sweep(clusters$centers, 2, col_sds, `*`)
  clusters$centers = sweep(clusters$centers, 2, col_means, `+`)
  clusterings[[k]] = clusters
}
plot(1:length(wss), wss, "l")

# plot per-cluster mean loss
mean_cluster_losses = aggregate(models$loss, by=list(clusterings[[6]]$cluster), FUN=mean)
for (i in 1:5) {
  cluster_losses = models[clusterings[[5]]$cluster == i, "loss"]
  plot(1:length(cluster_losses), cluster_losses, main=paste("Losses for cluster", i))
}