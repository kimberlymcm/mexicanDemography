# First sh file to run. Prepares the data


# Turns the vcf_ready file into the data_dict.
# Unclear where the vcf_ready file comes from atm
python2.7 ../../scripts/python/sfs_Alicia_Hernandez.py \
--vcf /srv/gsfs0/projects/bustamante/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_44M_vqsr_reheader_pass_20150826.vcf_ready \
--out /srv/gsfs0/projects/bustamante/kfm_projects/NA/NA_data/getIntrons/dadi_other/NA_CHB_exomes_20190301.data_dict \
--isField True

# Check that this new file is the same as the one I generated before in 2015. It is! Hooray.
diff /srv/gsfs0/projects/bustamante/kfm_projects/NA/NA_data/getIntrons/dadi_other/NA_CHB_exomes_20190301.data_dict /home/kfm/kfm_projects/NA/NA_data/getIntrons/remaskingNAH/NA_CHB_exomes_20150827.data_dict

