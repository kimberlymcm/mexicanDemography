# For bootstrapping part of the project
# Creates 1000 bootstrapped data_dict input files


reps <- 1000
blocks <- read.table("../../data/bootstrapping/bootstrap_blocks.txt")
colnames(blocks) <- c("chr", "start", "stop1", "SNP")
data <- read.table("/home/kfm/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_20150826.data_dict",
  header = TRUE)
data$chr <- gsub("chr", "", data$chr)
all_call <- vector()


# This could obviously be parallelized, but oh well for now.
for (i in 1:reps){

  print(i)
  # Random sample, with replacement
  new_blocks <- blocks[sample(nrow(blocks), nrow(blocks), replace=TRUE), ]
  result <- vector()
  callable_seq_len <- sum(new_blocks[, "SNP"])
  all_call <- c(all_call, callable_seq_len)
  
  # For each new block
  for (j in 1:NROW(new_blocks)){
    relevant_snps <- subset(data,
                     (chr == new_blocks[j, "chr"]) &
                     (pos <= new_blocks[j, "stop1"]) &
                     (pos > new_blocks[j, "start"])
                    )
    result <- rbind(result, relevant_snps)
  }

  result$num <- seq(1, NROW(result))
  fname <- paste(c("../../data/bootstrapping/input_bootstrap/Bootstrap_", i, ".txt"), collapse = "")
  write.table(result, fname, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

all_call = data.frame(all_call)
write.table(all_call, "../../data/bootstrapping/input_bootstrap/all_callable_lengths.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
