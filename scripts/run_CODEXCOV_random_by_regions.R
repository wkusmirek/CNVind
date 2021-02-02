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

  indices <- sample(1:nrow(Y), number_of_neighbours, replace=F)
  indices <- sort(unique(c(target_id,indices)))

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
  write.table(t(Yhat_1[which(indices==target_id),]), paste(output_Yhat_file,'_1',sep=''), sep=',', row.names=FALSE, col.names=FALSE)
  write.table(t(Yhat_2[which(indices==target_id),]), paste(output_Yhat_file,'_2',sep=''), sep=',', row.names=FALSE, col.names=FALSE)
  write.table(t(Yhat_3[which(indices==target_id),]), paste(output_Yhat_file,'_3',sep=''), sep=',', row.names=FALSE, col.names=FALSE)
}
