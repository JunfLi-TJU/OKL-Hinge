setwd("D:/experiment/AAAI2023/code")
rm(list = ls())
library(MASS)

d_index <- 2

dpath          <- file.path("D:/experiment/AAAI2023/dataset/")  

Dataset        <- c("w8a", "magic04-N", "ijcnn1_all","a9a_all","SUSY50000","mushrooms","mnist12")

savepath      <- paste0("D:/experiment/AAAI2023/Result/",paste0("FOGD-",Dataset[d_index],".txt"))

traindatapath    <- file.path(dpath, paste0(Dataset[d_index], ".train"))                
traindatamatrix  <- as.matrix(read.table(traindatapath))
trdata           <- traindatamatrix[ ,-1]
ylabel           <- traindatamatrix[ ,1]

length_tr        <- nrow(trdata)    
feature_tr       <- ncol(trdata)  

U <- 15
p_setting <-list(
  gamma = 4,
  d     = feature_tr,
  D     = 400,
  zx    = array(0,dim=2*400),
  eta   = 100/sqrt(length_tr)
) 

reptimes <- 10

runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
cum_Loss  <- c(rep(0, reptimes))
tem_data  <- c(rep(0,feature_tr))
norm_w     <- c(rep(0, reptimes))

for(re in 1:reptimes)
{
  order <- sample(1:length_tr,length_tr,replace = F)   #dis
  error <- 0
  w     <- c(rep(0, 2*p_setting$D))
  sigma <- (1/p_setting$gamma)^2 *diag(p_setting$d)
  u     <- mvrnorm(p_setting$D,rep(0,p_setting$d),sigma)   # w--->D*d
  t1    <- proc.time()                                     #proc.time()
  
  for (i in 1:length_tr)
  {
    tem = u%*%trdata[order[i],]
    coszx <- cos(tem)/sqrt(p_setting$D)
    sinzx <- sin(tem)/sqrt(p_setting$D)
    p_setting$zx <- c(coszx,sinzx)
    sum   <- crossprod(w,p_setting$zx)[1,1]
    hat_y <- 1
    if(sum < 0)
      hat_y <- -1
    if(hat_y != ylabel[order[i]])
      error = error+1
    if(ylabel[order[i]]*sum<1)
    {
      w     <- w + p_setting$eta*ylabel[order[i]]*p_setting$zx 
    }
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  norm_w[re] <- sqrt(crossprod(w,w)[1,1])
}

save_result <- list(
  note     = c(" the next term are:alg_name--dataname--ker_para--sam_num--RFF--run_time--tot_run_time--ave_run_time--err_num--all_err_rate--ave_err_rate--sd_time--sd_err"),
  alg_name = c("FOGD"),
  dataname = paste0(Dataset[d_index], ".train"),
  ker_para = p_setting$gamma,
  sam_num  = length_tr,
  rff_num  = p_setting$D,
  run_time = as.character(runtime),
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  err_num  = errorrate,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)
write.table(save_result,file=savepath,row.names =TRUE, col.names =FALSE, quote = T)

sprintf("the kernel parameter is %f", p_setting$gamma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of random features are %d", p_setting$D)
sprintf("total running time is %.1f in dataset", sum(runtime))
sprintf("average running time is %.3f in dataset", sum(runtime)/reptimes)
sprintf("the average MR is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of MR is %.5f in dataset", sd(errorrate))
sprintf("average norm is %.5f in dataset", mean(norm_w))
