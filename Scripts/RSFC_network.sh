# USE THIS TO EXTRACT THE PE VALUES FOR EACH NETWORK AND EACH CONDITION

cohort="INSIGHT"

dual_type="Yeo7NStriatbinary_17_1_DMTvsPCBpostNoRat_GSR/"
thresh="thr5"

cohort_dir="/Users/christophertimmermann/Documents/Imaging_fMRI/"
DUALdir="$cohort_dir/Raw/analysis/$dual_type"
DUALcorr_dir="$DUALdir/withinRSN" ; mkdir -p $DUALcorr_dir

# SELECT YEO 7 NETWORKS
ICA_good="00 01 02 03 04 05 06"



# GENERATE A FOLDER WHICH OUTPUTS THE MEAN PE VALUES PER NETWORK PER SUBJECT PER CONDITION< VALUES ARE ORDERED IN THE SAME WAY AS THE INPUT FILES FOR DUAL_REGRESSION.SH SCRIPT
touch ${DUALcorr_dir}/Intra_RSN${RSN}.txt

do_fixed=1
if [[ "$do_fixed" == 1 ]]; then  #for 15 subjectes
	for ICA in `echo $ICA_good`
	do
		mkdir -p $DUALdir/Fixed_ICA${ICA}
		cp $DUALdir/dr_stage2_ic00${ICA}.nii.gz $DUALdir/Fixed_ICA${ICA}/dr_stage2_ic00${ICA}.nii.gz
		pushd $DUALdir/Fixed_ICA${ICA}; fslsplit dr_stage2_ic00${ICA}.nii.gz; popd
	done
fi


do_INTRA_mid_and_After=1
if [[ "$do_INTRA_mid_and_After" == 1 ]]; then
	for ICA in `echo $ICA_good`
	do
		for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 #One for each file


		do
			touch ${DUALcorr_dir}/INTRA_mid_and_After_${ICA}.txt
			fslmeants -i $DUALdir/Fixed_ICA${ICA}/vol00${i}.nii.gz -m /Users/christophertimmermann/Documents/Imaging_fMRI/Scripts/New/Dual_reg/Yeo7NStriatBinaryNets/vol00${ICA}.nii.gz -o ${DUALcorr_dir}/out.txt

			cat ${DUALcorr_dir}/out.txt >> ${DUALcorr_dir}/INTRA_mid_and_After_${ICA}.txt
			rm ${DUALcorr_dir}/out.txt
		done
	done
fi

















