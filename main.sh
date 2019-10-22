#!/bin/bash

 

export OMP_NUM_THREADS=16

export MKL_NUM_THREADS=16

export DYLD_FALLBACK_LIBRARY_PATH=/imaging/mlr_imaging/AH/afni/v18.3.03/

export AFNI_3dDespike_NEW=YES

FSLOUTPUTTYPE=NIFTI_GZ

 

#set generic variables manually

ROOT=~/Desktop/mlr-lab/AH/Projects/AHalai/Ongoing/P1370/

#FOLDERS="$ROOT"/*

#Label of final output file, typically name of experiment

exp=SemanticLoc

TR=1500

fwhm=4

cd "$ROOT"

 

for subj in "$ids"; do

 

for cond in Run1 Run2 Run3; do

cd "$ROOT""$subj"

 

# Check if T1.anat folder exists, if not run fsl_anat to preprocess T1

echo +* "run fsl_anat on all T1 datasets"

if [ -d "T1.anat" ]; then

echo +* "T1.anat folder already exists - not running fsl_anat"

else

fsl_anat -i t1_mprage_sag_1 -o T1 --nosubcortseg -t T1

fi

 

echo +* "Perform pre-processing $cond"

filename="$exp""$cond"

 

mkdir -p "$filename"

cp "$ROOT"timing.txt "$filename"

cd "$filename"/

 

echo "$ROOT""$subj"_"$filename"

 

echo "

++++++++++++++++++++++++"

echo +* "Set up script run environment"

echo "

++++++++++++++++++++++++"

echo +* "Copy in functional datasets, reset NIFTI tags as needed"

3dcalc -a ../"$filename"_1.nii.gz -expr 'a' -prefix ./"$filename"_1.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./"$filename"_1.nii -overwrite

3dcalc -a ../"$filename"_2.nii.gz -expr 'a' -prefix ./"$filename"_2.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./"$filename"_2.nii -overwrite

3dcalc -a ../"$filename"_3.nii.gz -expr 'a' -prefix ./"$filename"_3.nii

nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./"$filename"_3.nii -overwrite

 

echo "

++++++++++++++++++++++++"

echo +* "Calculate and save motion parameters"

if [ $cond = "Run1" ];then

3dDespike -overwrite -prefix ./"$filename"_1_vrA.nii.gz ./"$filename"_1.nii

3daxialize -overwrite -prefix ./"$filename"_1_vrA.nii.gz ./"$filename"_1_vrA.nii.gz

3dcalc -a ./"$filename"_1_vrA.nii.gz[0]  -expr 'a' -prefix "$exp"base.nii.gz

3dvolreg -overwrite -prefix ./"$filename"_1_vrA.nii.gz -base "$exp"base.nii.gz -dfile ./"$filename"_1_vrA.1D -1Dmatrix_save ./"$filename"_1_vrmat.aff12.1D ./"$filename"_1_vrA.nii.gz

1dcat ./"$filename"_1_vrA.1D[1..6] > "$filename"_motion.1D

cat "$filename"_motion.1D > "$filename"_motion.txt

#copy base image to root folder as this will be used as the image to align all runs

cp ./"$exp"base.nii.gz ../"$exp"base.nii.gz

else

3dDespike -overwrite -prefix ./"$filename"_1_vrA.nii.gz ./"$filename"_1.nii

3daxialize -overwrite -prefix ./"$filename"_1_vrA.nii.gz ./"$filename"_1_vrA.nii.gz

3dvolreg -overwrite -prefix ./"$filename"_1_vrA.nii.gz -base ../"$exp"base.nii.gz -dfile ./"$filename"_1_vrA.1D -1Dmatrix_save ./"$filename"_1_vrmat.aff12.1D ./"$filename"_1_vrA.nii.gz

1dcat ./"$filename"_1_vrA.1D[1..6] > "$filename"_motion.1D

cat "$filename"_motion.1D > "$filename"_motion.txt

fi

 

echo "

++++++++++++++++++++++++"

echo +* "Preliminary preprocessing of functional datasets: despike, tshift, deoblique, and/or axialize"

echo --------"Preliminary preprocessing dataset "$filename"_1.nii.gz of TE=12ms to produce e1_ts+orig"

3dDespike -overwrite -prefix ./"$filename"_1_pt.nii.gz "$filename"_1.nii

3dTshift -heptic -tpattern @timing.txt -prefix ./e1_ts+orig ./"$filename"_1_pt.nii.gz

3drefit -view orig e1_ts*HEAD

3daxialize  -overwrite -prefix ./e1_ts+orig ./e1_ts+orig

3drefit -deoblique -TR $TR e1_ts+orig

echo --------"Preliminary preprocessing dataset "$filename"_2.nii.gz of TE=24.83ms to produce e2_ts+orig"

3dDespike -overwrite -prefix ./"$filename"_2_pt.nii.gz "$filename"_2.nii

3dTshift -heptic -tpattern @timing.txt -prefix ./e2_ts+orig ./"$filename"_2_pt.nii.gz

3drefit -view orig e2_ts*HEAD

3daxialize  -overwrite -prefix ./e2_ts+orig ./e2_ts+orig

3drefit -deoblique -TR $TR e2_ts+orig

echo --------"Preliminary preprocessing dataset "$filename"_3.nii.gz of TE=37.66ms to produce e3_ts+orig"

3dDespike -overwrite -prefix ./"$filename"_3_pt.nii.gz "$filename"_3.nii

3dTshift -heptic -tpattern @timing.txt -prefix ./e3_ts+orig ./"$filename"_3_pt.nii.gz

3drefit -view orig e3_ts*HEAD

3daxialize  -overwrite -prefix ./e3_ts+orig ./e3_ts+orig

3drefit -deoblique -TR $TR e3_ts+orig

 

if [ $cond = "Run1" ]; then

echo "

++++++++++++++++++++++++"

echo +* "Preparing functional masking for this ME-EPI run."

fslmaths "$filename"_1 -Tmean mean

flirt -in mean -ref ../T1.anat/T1_biascorr_brain -omat fun2str.mat -cost normmi -dof 6 -out cr_mean_T1 -interp nearestneighbour

convert_xfm -omat str2fun.mat -inverse fun2str.mat

flirt -in ../T1.anat/T1_biascorr_brain -ref mean -applyxfm -init str2fun.mat -out mask

fslmaths mask -fillh -bin -nan mask

 

cp ./mask.nii.gz ../"$exp"mask.nii.gz

fi

 

cp ../"$exp"mask.nii.gz ./mask.nii.gz

 

echo --------"Apply combined normalization/co-registration/motion correction parameter set to e1_ts+orig"

3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply "$filename"_1_vrmat.aff12.1D -base mask.nii.gz -input e1_ts+orig -prefix ./e1_vr.nii.gz

3dcalc -float -overwrite -a mask.nii.gz -b ./e1_vr.nii.gz -expr 'step(a)*b' -prefix ./e1_in.nii.gz

3dTstat -prefix ./e1_mean.nii.gz ./e1_in.nii.gz

3dTstat -stdev -prefix ./e1_std.nii.gz ./e1_in.nii.gz

echo --------"Apply combined normalization/co-registration/motion correction parameter set to e2_ts+orig"

3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply "$filename"_1_vrmat.aff12.1D -base mask.nii.gz -input e2_ts+orig -prefix ./e2_vr.nii.gz

3dcalc -float -overwrite -a mask.nii.gz -b ./e2_vr.nii.gz -expr 'step(a)*b' -prefix ./e2_in.nii.gz

3dTstat -prefix ./e2_mean.nii.gz ./e2_in.nii.gz

3dTstat -stdev -prefix ./e2_std.nii.gz ./e2_in.nii.gz

echo --------"Apply combined normalization/co-registration/motion correction parameter set to e3_ts+orig"

3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply "$filename"_1_vrmat.aff12.1D -base mask.nii.gz -input e3_ts+orig -prefix ./e3_vr.nii.gz

3dcalc -float -overwrite -a mask.nii.gz -b ./e3_vr.nii.gz -expr 'step(a)*b' -prefix ./e3_in.nii.gz

3dTstat -prefix ./e3_mean.nii.gz ./e3_in.nii.gz

3dTstat -stdev -prefix ./e3_std.nii.gz ./e3_in.nii.gz

 

source activate AHALAI

tedana -d e1_in.nii.gz e2_in.nii.gz e3_in.nii.gz -e 12 24.83 37.66 --png --mask mask.nii.gz --maxit 500 --maxrestart 10

source deactivate

 

if [ $cond = "Run1" ]; then

echo "

++++++++++++++++++++++++"

echo +* "Coregister tedana outputs to T1 using FLIRT."

fslmaths e1_in -Tmean -thr 0 -nan mean

fslreorient2std mean mean

flirt -in mean -ref ../T1.anat/T1_biascorr_brain -omat fun2str.mat -cost normmi -dof 6 -out cr_mean_T1 -interp nearestneighbour

convert_xfm -omat str2fun.mat -inverse fun2str.mat

fi

 

echo +* "Change datatype to float to reduce file size"

gzip -f ts_OC.nii dn_ts_OC.nii r_ts_OC.nii r_dn_ts_OC.nii betas_OC.nii lowk_ts_OC.nii hik_ts_OC.nii mepca_OC_components.nii

 

fslmaths ts_OC -mul 1 ts_OC -odt float

fslmaths dn_ts_OC -mul 1 dn_ts_OC -odt float

fslmaths r_ts_OC -mul 1 r_ts_OC -odt float

fslmaths r_dn_ts_OC -mul 1 r_dn_ts_OC -odt float

fslmaths betas_OC -mul 1 betas_OC -odt float

fslmaths lowk_ts_OC -mul 1 lowk_ts_OC -odt float

fslmaths hik_ts_OC -mul 1 hik_ts_OC -odt float

fslmaths mepca_OC_components -mul 1 mepca_OC_components -odt float

 

echo +* "Normalising medn and tsoc to MNI space using FNIRT warps"

fslreorient2std dn_ts_OC r_dn_ts_OC

fslreorient2std ts_OC r_ts_OC

 

applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm -i r_dn_ts_OC -w ../T1.anat/T1_to_MNI_nonlin_coeff --premat=../"$exp"Run1/fun2str.mat -o wcr_dn_ts_OC --interp=spline

applywarp -r $FSLDIR/data/standard/MNI152_T1_2mm -i r_ts_OC -w ../T1.anat/T1_to_MNI_nonlin_coeff --premat=../"$exp"Run1/fun2str.mat -o wcr_ts_OC --interp=spline

 

fslmaths wcr_dn_ts_OC -mul $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil1 -thr 0 -nan wcr_dn_ts_OC

fslmaths wcr_ts_OC -mul $FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil1 -thr 0 -nan wcr_ts_OC

 

echo +* "smooth"

sigma=$(echo "$fwhm / 2.35482004503" | bc -l)

echo $sigma

 

fslmaths wcr_dn_ts_OC -kernel gauss $sigma -fmean -nan swcr_dn_ts_OC

fslmaths wcr_ts_OC -kernel gauss $sigma -fmean -nan swcr_ts_OC

 

done

 

done