# Calculate 500 kb blocks based on the bed files
# Input:
#    chr_lens: Length of each chromosome
#    targets: Locations of target sites in bed file format
# Output: 
#    blocks: A file with each row as
#            chr start_pos_of_block end_pos_of_block #_sites_in_block


block_size <- 500000
chr_lens <- read.table("/home/kfm/kfm_projects/NA/NA_data/getIntrons/dadi_other/chr_lengths.txt")
#targets <- read.table('/home/kfm/kfm_projects/NA/NA_data/compareCaptureRegions/TargetRegions_Agilent_44M_v2_reformat.list.bed.noChr') #chr, start (zero-based), stop (excluded) # didn't remove this in case knowing its location is useful for later
targets <- read.table("/home/kfm/kfm_projects/NA/NA_data/getCallableLength/4fold_intron_outgroup_20150826.sites.only.bed")
colnames(targets) <- c("chr", "start", "pos")
targets$chr <- gsub("chr", "", targets$chr)
blocks <- vector()
total_blocks <- 0


for (j in 1:22) {
  print(j)
  targ_chr <- subset(targets, chr == j) # Limit to sites in relevant chr
  len <- chr_lens[j, 2]
  num_blocks <- round(len / block_size) # Num blocks to break chr into
  start_pos <- 0
  for (i in 1:num_blocks){
    block_num <- i
    end_pos <- start_pos + block_size
    targ_in_subset <- subset(targ_chr, pos <= end_pos & pos > start_pos)
    if (nrow(targ_in_subset) != 0){
      blocks <- rbind(blocks, c(j, start_pos, end_pos, nrow(targ_in_subset)))
    }
    start_pos <- end_pos
    total_blocks <- num_blocks + total_blocks
  }
}
print(total_blocks)
write.table(blocks, "../../data/bootstrapping/bootstrap_blocks.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
