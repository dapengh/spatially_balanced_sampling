source(file = "./functions.R")
packages <- c("spsurvey", "BalancedSampling", "parallel", "foreach", "doParallel", "doRNG")
lapply(packages, library, character.only = TRUE)


t=25
seed = 202117455
cl <- makeCluster(16)
registerDoParallel(cl)
set.seed(20200326)
tmp <- foreach(i = 1:100, .errorhandling = "remove") %dorng% {
  spatial_sampling(time = t, initial_prev = 0.1, rate = 0.2, sample_size = 60, initial_loc = "m", cut_off_prev = 0.1, seed = seed,n_iter = 100)
}
stopCluster(cl)

pd <- data.frame(Reduce("+", tmp) / length(tmp))
names(pd) <-  c("cube", "lcube", "lpm", "scps", "grts", "srs", "prev")
write.csv(pd, file = "./res/pd_ini010_rate02_size60_m.csv", row.names = F)