""" Converts from some format 'vcf_ready' to dadi input format """

import os
import argparse

import numpy
from numpy import array

import dadi
#import pylab
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec

USAGE = """
sfs.py      --vcf <vcf>
            --data_dict <data file in dadi format>
            --isField <determines whether ancestry calls are a field in the VCF INFO field>
            --out <output file>
"""
parser = argparse.ArgumentParser(USAGE)

parser.add_argument("--vcf", default="/srv/gsfs0/projects/bustamante/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_44M_vqsr_reheader_pass_20150826.vcf_ready")
parser.add_argument("--out", default="/srv/gsfs0/projects/bustamante/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_20150826.data_dict")
parser.add_argument("--data", default="/srv/gsfs0/projects/bustamante/kfm_projects/NA/dadi/data/NatAmExomes.gape.syn.mask.auto.data_dict")
parser.add_argument("--isField", default=True, type=bool)
args = parser.parse_args()

def count_ref_alt(snpLine, ind_indices):
    count_ref = 0
    count_alt = 0
    for ind in ind_indices:
        #print ind
        if snpLine[ind].startswith("0/0"): #check this one
            count_ref = count_ref + 2
        elif snpLine[ind].startswith("0/1") or snpLine[ind].startswith("1/0"):
            count_ref = count_ref + 1
            count_alt = count_alt + 1
        elif snpLine[ind].startswith("1/1"):
            count_alt = count_alt + 2
    return(str(count_ref), str(count_alt))


if (args.data.startswith("/home") or args.data.startswith("/srv")): # we need to create the fs file
    vcf = open(args.vcf)
    out = open(args.out, 'w')

    out.write("Human\tChimp\tAllele1\tHUI\tMYA\tNAH\tTAR\tTRQ\tCHB\tAllele2\tHUI\tMYA\tNAH\tTAR\tTRQ\tCHB\tchr\tpos\n")
    for line in vcf:
        if line.startswith("##"):
            pass
        elif line.startswith("#"): # This is the one line with all the stuff on it
            header = line.strip().split()
            tar = []
            hui = [] 
            mya = []
            trq = []
            nah = []
            chb = []
            tar_h = []
            hui_h = []
            mya_h = []
            trq_h = []
            nah_h = []
            chb_h = []
            for ind in header[9:len(header)]:
                if ind.startswith("TAR"):
                    tar.append(header.index(ind))
                    tar_h.append(header[header.index(ind)])
                elif ind.startswith("HUI"):
                    hui.append(header.index(ind))
                    hui_h.append(header[header.index(ind)])
                elif ind.startswith("MYA"):
                    mya.append(header.index(ind))
                    mya_h.append(header[header.index(ind)])
                elif ind.startswith("TRQ"):
                    trq.append(header.index(ind))
                    trq_h.append(header[header.index(ind)])
                elif ind.startswith("NA18"):
                    chb.append(header.index(ind))
                    chb_h.append(header[header.index(ind)])
                else: # It is NAH
                    nah.append(header.index(ind))
                    nah_h.append(header[header.index(ind)])
            print 'TAR: ' + str(len(tar))
            print tar
            print tar_h
            print 'HUI: ' + str(len(hui))
            print hui
            print hui_h
            print 'MYA: ' + str(len(mya))
            print mya
            print mya_h
            print 'TRQ: ' + str(len(trq))
            print trq
            print trq_h
            print 'NAH: ' + str(len(nah))
            print nah
            print nah_h
            print 'CHB: ' + str(len(chb))
            print chb
            print chb_h
        else:
            myLine = line.strip().split()
            if args.isField: # This does make a bunch of sense. isField=1, so I guess it always exists?
                anc_field = str(myLine[7].split("Tri_hg19=")[1].split(";")[0].upper().strip()) # ;Tri_hg19=CCG;
                anc_field2 = str(myLine[7].split(";Tri_panTro4=")[1].split(";")[0].upper().strip()) # ;Tri_hg19=CCG;

                if myLine[4] == "A" or myLine[4] == "C" or myLine[4] == "G" or myLine[4] == "T": # This means not triallelic
                    if anc_field[1] == myLine[3] or anc_field[1] == myLine[4]: # This means that one of the alleles genotyped is the reference
                        if anc_field[1] == myLine[3]: # Idk if this really matters, but it makes sure allele1 is the correct one
                            first = myLine[3]
                            second = myLine[4]
                            (count_ref_hui, count_alt_hui) = count_ref_alt(myLine, hui)
                            (count_ref_mya, count_alt_mya) = count_ref_alt(myLine, mya)
                            (count_ref_nah, count_alt_nah) = count_ref_alt(myLine, nah)
                            (count_ref_tar, count_alt_tar) = count_ref_alt(myLine, tar)
                            (count_ref_trq, count_alt_trq) = count_ref_alt(myLine, trq)
                            (count_ref_chb, count_alt_chb) = count_ref_alt(myLine, chb)
                        else:
                            first = myLine[4]
                            second = myLine[3]
                            (count_alt_hui, count_ref_hui) = count_ref_alt(myLine, hui)
                            (count_alt_mya, count_ref_mya) = count_ref_alt(myLine, mya)
                            (count_alt_nah, count_ref_nah) = count_ref_alt(myLine, nah)
                            (count_alt_tar, count_ref_tar) = count_ref_alt(myLine, tar)
                            (count_alt_trq, count_ref_trq) = count_ref_alt(myLine, trq)
                            (count_alt_chb, count_ref_chb) = count_ref_alt(myLine, chb)
                        out.write(anc_field + "\t" + anc_field2 + "\t" + first + "\t" + count_ref_hui + "\t" + count_ref_mya + 
                            "\t" + count_ref_nah + "\t" + count_ref_tar + "\t" + count_ref_trq + "\t" + count_ref_chb + "\t" +
                            second + "\t" + count_alt_hui + "\t" + count_alt_mya + "\t" + count_alt_nah + "\t" + count_alt_tar +
                            "\t" + count_alt_trq + "\t" + count_alt_chb + "\t" + myLine[0] + "\t" + myLine[1] + "\n")

                    else:
                        pass
                else:
                    pass

    
else: #just load the data dict file
    data_dict = dadi.Misc.make_data_dict(options.data)
    hui_mya = dadi.Spectrum.from_data_dict(data_dict, ['HUI', 'MYA'], [26, 26], polarized=True)
    hui_nah = dadi.Spectrum.from_data_dict(data_dict, ['HUI', 'NAH'], [26, 32], polarized=True)
    hui_tar = dadi.Spectrum.from_data_dict(data_dict, ['HUI', 'TAR'], [26, 38], polarized=True)
    hui_trq = dadi.Spectrum.from_data_dict(data_dict, ['HUI', 'TRQ'], [26, 30], polarized=True)
    mya_nah = dadi.Spectrum.from_data_dict(data_dict, ['MYA', 'NAH'], [26, 32], polarized=True)
    mya_tar = dadi.Spectrum.from_data_dict(data_dict, ['MYA', 'TAR'], [26, 38], polarized=True)
    mya_trq = dadi.Spectrum.from_data_dict(data_dict, ['MYA', 'TRQ'], [26, 30], polarized=True)
    nah_tar = dadi.Spectrum.from_data_dict(data_dict, ['NAH', 'TAR'], [32, 38], polarized=True)
    nah_trq = dadi.Spectrum.from_data_dict(data_dict, ['NAH', 'TRQ'], [32, 30], polarized=True)
    tar_trq = dadi.Spectrum.from_data_dict(data_dict, ['TAR', 'TRQ'], [38, 30], polarized=True)
    chb_hui = dadi.Spectrum.from_data_dict(data_dict, ['CHB', 'HUI'], [103, 26], polarized=True)
    chb_mya = dadi.Spectrum.from_data_dict(data_dict, ['CHB', 'MYA'], [103, 26], polarized=True)
    chb_nah = dadi.Spectrum.from_data_dict(data_dict, ['CHB', 'NAH'], [103, 32], polarized=True)
    chb_tar = dadi.Spectrum.from_data_dict(data_dict, ['CHB', 'TAR'], [103, 38], polarized=True)
    chb_trq = dadi.Spectrum.from_data_dict(data_dict, ['CHB', 'TRQ'], [103, 30], polarized=True)
    
    
    #subplots_adjust(wspace=0.2,hspace=0.2)
    #
    #gs = gridspec.GridSpec(5, 2, wspace=0.5, hspace=0.3)
    #plt.subplot(gs[0, 0])
    fig = plt.figure()
    fig.set_size_inches(6, 10)
    ax1 = fig.add_subplot(5,2,1, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(hui_mya, vmin=0.1)
    #plt.subplot(gs[0, 1])
    ax2 = fig.add_subplot(5,2,2, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(hui_nah, vmin=0.1)
    ax3 = fig.add_subplot(5,2,3, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(hui_tar, vmin=0.1)
    ax4 = fig.add_subplot(5,2,4, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(hui_trq, vmin=0.1)
    ax5 = fig.add_subplot(5,2,5, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(mya_nah, vmin=0.1)
    ax6 = fig.add_subplot(5,2,6, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(mya_tar, vmin=0.1)
    ax7 = fig.add_subplot(5,2,7, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(mya_trq, vmin=0.1)
    ax8 = fig.add_subplot(5,2,8, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(nah_tar, vmin=0.1)
    ax9 = fig.add_subplot(5,2,9, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(nah_trq, vmin=0.1)
    ax10 = fig.add_subplot(5,2,10, aspect=1)
    dadi.Plotting.plot_single_2d_sfs(tar_trq, vmin=0.1)
    pylab.show()
    pylab.savefig(options.out + '_unfolded.pdf')

