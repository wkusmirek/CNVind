run_CODEXCOV <- function(K_from,
                         K_to,
                         lmax,
                         number_of_neighbours,
                         target_id,                         
                         input_cov_table,
                         input_bed,
                         output_Yhat_file){

  Y <- read.csv(input_cov_table)
  targets <- read.delim(input_bed)
  rownames(Y) <- 1:nrow(Y)
  rownames(targets) <- 1:nrow(targets)

  correlation <- cor(t(Y[target_id,]),t(Y[,]))
  thr <- sort(correlation,decreasing=TRUE)[number_of_neighbours]
  indices <- which(correlation >= thr)
  indices <- as.numeric(gsub("\\D+", "", indices))

  print(indices)

  Y <- Y[indices,]
  targets <- targets[indices,]  
  rownames(Y) <- 1:nrow(Y)
  rownames(targets) <- 1:nrow(targets)

  ref <- IRanges(start = targets[,"st_bp"], end = targets[,"ed_bp"])
  gcmapp1_result <- gcmapp1(targets[1,'chr'], ref)
  gc <- gcmapp1_result$gc

  normObj_result <- normObj1(as.matrix(Y), gc, K=K_from:K_to)
  Yhat <- normObj_result$Yhat
  AIC <- normObj_result$AIC
  BIC <- normObj_result$BIC
  RSS <- normObj_result$RSS
  K <- normObj_result$K
  optK <- K[which.max(BIC)]
  Yhat_opt <- Yhat[[which(K == optK)]][,]
  Yhat_1 <- Yhat[[which(K == 1)]][,]
  Yhat_2 <- Yhat[[which(K == 2)]][,]
  Yhat_3 <- Yhat[[which(K == 3)]][,]
  write.table(t(Yhat_opt[which(indices==target_id),]), output_Yhat_file, sep=',', row.names=FALSE, col.names=FALSE)
}
