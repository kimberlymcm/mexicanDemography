# Code for calculating confidence intervals for bootstrapped results
# Copied from google doc:
#   Dadi_Results/Bootstrapping/Procedure doc

prev <- 1
for (i in 1:1000){

    fname <- paste(c("../../results/bootstrapping/results_0.09Tb/20190224_bootstrap.", i),
        collapse = "")

    # In the google doc, the below line was incorrectly placed here.
    # So I moved it below. Luckily, this way (moved below) correctly replicates the results
    # The order that we read in the bootstrapped results matters because we match it up with the call length  
    # if (!exists("dataset")){ dataset <- read.table(fname) }
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
          temp_dataset <- read.table(fname)
          dataset <- rbind(dataset, temp_dataset[1, ])
          rm(temp_dataset)
        if((nrow(dataset) - prev ) != 1){
            print(nrow(dataset) - prev)
            print(nrow(dataset))
        }
        prev <- nrow(dataset)
    }
    if (!exists("dataset")){ dataset <- read.table(fname) }
}

a <- read.table("../../data/bootstrapping/input_bootstrap/all_callable_lengths.txt")

N_anc <- dataset[, 2] / (4 * a * 0.0000000125 * 0.96)
quantile(N_anc$V1, c(0.025, 0.5, 0.975))
mean(N_anc$V1)

Nb <- N_anc$V1 * dataset[, 4]
mean(Nb)
quantile(Nb, c(0.025, 0.5, 0.975))

N_Tar <- N_anc$V1 * dataset[, 8]
mean(N_Tar)
quantile(N_Tar, c(0.025, 0.5, 0.975))

N_Hui <- N_anc$V1 * dataset[, 9]
mean(N_Hui)
quantile(N_Hui, c(0.025, 0.5, 0.975))

N_Trq <- N_anc$V1 * dataset[, 10]
mean(N_Trq)
quantile(N_Trq, c(0.025, 0.5, 0.975))

N_Mya <- N_anc$V1 * dataset[, 11]
mean(N_Mya)
quantile(N_Mya, c(0.025, 0.5, 0.975))

T_TrqMya <- N_anc$V1 * dataset[, 7] * 2 * 29 
mean(T_TrqMya)
quantile(T_TrqMya, c(0.025, 0.5, 0.975))

T_TarHui <- N_anc$V1 * dataset[, 6] * 2 * 29 + T_TrqMya
mean(T_TarHui)
quantile(T_TarHui, c(0.025, 0.5, 0.975))

T_s1 <- N_anc$V1 * dataset[, 5] * 2 * 29 + T_TarHui 
mean(T_s1)
quantile(T_s1, c(0.025, 0.5, 0.975))

Tb <- N_anc$V1 * dataset[,3] * 2 * 29 + T_s1
mean(Tb)
quantile(Tb, c(0.025, 0.5, 0.975))
