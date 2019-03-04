"""
   NOTE: Called 'mask_eur_Alicia_EDITED.py' elsewhere 
"""

import argparse
import gzip


# TODO(kmcmanus): Turn this into a main method
parser = argparse.ArgumentParser(description="Mask european ancestry in vcfs")

parser.add_argument("--vcf", default="/srv/gs1/projects/bustamante/bgi_data/Native_American_exomes/processed/NA_CHB_exomes_44M_vqsr_reheader.vcf.gz")
parser.add_argument("--bed_files", default="/srv/gs1/projects/bustamante/armartin_projects/NA_analysis/nat_mask/bed_list_20150822.txt")
parser.add_argument("--out", default="/srv/gs1/projects/bustamante/bgi_data/Native_American_exomes/processed/NA_CHB_exomes_44M_vqsr_reheader_NOTpass_mask_20150902.vcf")

args = parser.parse_args()

if args.vcf.endswith("gz"):
    vcf = gzip.open(args.vcf)
    #out = gzip.open(options.out, "w")
    out = open(args.out, "w")
else:
    vcf = open(args.vcf)
    out = open(args.out, "w")

bed_files = open(args.bed_files)
# filter_inds appears to just include all individuals?
filter_inds = set()
for line in bed_files:
    filter_inds.add(line.strip().split(".")[0])

chrs = map(str, range(1, 23))

#requires a chr, saves all position starts and ends to mask for all individuals
def chr_pos(chr):
    bed_files = open(args.bed_files)
    ind_starts = {}
    ind_ends = {}
    for line in bed_files:
        ind_bed = line.strip()
        #print ind_bed
        ind = ind_bed.split(".")[0]
        ind_mask = open("/srv/gs1/projects/bustamante/armartin_projects/NA_analysis/nat_mask/" + ind_bed)
        for mask in ind_mask:
            mask_line = mask.strip().split()
            if mask_line[0] == chr:
                if ind in ind_starts:
                    ind_starts[ind].append(mask_line[1])
                    ind_ends[ind].append(mask_line[2])
                else:
                    ind_starts[ind] = [mask_line[1]]
                    ind_ends[ind] = [mask_line[2]]
    return(ind_starts, ind_ends)

# I don't think the below line is used for anything?
# (ind_starts, ind_ends) = chr_pos("1")

#determine whether position is between start and end intervals for that ind (bed files are NAT_NAT regions)
def ind_pos(position, ind, current_geno, chr_starts, chr_ends):
    """ Masks the genotype call if it is not in a native segment.
        It does this by determining whether position is between start and end
        intervals for that ind (bed files are NAT_NAT regions
    """
    ind_starts = chr_starts[ind]
    ind_ends = chr_ends[ind]
    #print [position, ind, current_geno, ind_starts, ind_ends]
    in_interval = False
    for interval in range(len(ind_starts)):
        if position > int(ind_starts[interval]) and position < int(ind_ends[interval]):
            in_interval = True
            break
    if in_interval:
        return(current_geno)
    else:
        return("./.")


last_chr = "0"
for line in vcf:
    myLine = line.strip()
    if myLine.startswith("#"):
        out.write(myLine + "\n")
        if myLine.startswith("#CHROM"):
            header = myLine.split()
    else:
        myLine = myLine.split()
        current_chr = myLine[0].split("chr")[1]
        # Presumbly we are ignoring the X and Y
        if current_chr != "X" and current_chr != "Y":
            if current_chr != last_chr:
                (chr_starts, chr_ends) = chr_pos(current_chr)
            out.write("\t".join(myLine[0:9]) + "\t")
            current_pos = int(myLine[1])
            to_write = []
            for ind in header[9:len(header)]:
                # In this case, I think all inds need to be filtered. However, this if/else statement is just in
                # case there are inds in the vcf that aren't listed in the bed files
                if ind in filter_inds:
                    to_write.append(ind_pos(current_pos, ind, myLine[header.index(ind)], chr_starts, chr_ends))
                else:
                    print myLine[header.index(ind)]
                    to_write.append(myLine[header.index(ind)])
            out.write("\t".join(to_write))
            out.write("\n")
            last_chr = current_chr

