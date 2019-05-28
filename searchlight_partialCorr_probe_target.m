function searchlight_partialCorr_probe_target
clear;
clc;




basedir = '/scratch/groups/Projects/P1369/zygao/fMRI_data';
%subs=1:20;

subs=1;% for debug


epsilon =  1e-6;

xlength =  91;
ylength =  109;
zlength =  91;
radius  =  2;     % the cubic size is (2*radius+1) by (2*radius+1) by (2*radius+1)
step    =  1;

tic;


load(sprintf('%s/beh_results/sem/allw2b.mat',basedir));



for iSub=1:length(subs)
    sub_current=subs(iSub);
    fprintf('Subject %d is processing.\n',sub_current);
    
    load(sprintf('%s/beh_results/sem/sub%02d_alldata.mat',basedir,sub_current));
    
    ifg=[];pmtg=[];rpmtg=[];
    for run=1:4
        ifg_tmp=load(sprintf('%/mps/stand/roiTP/sub%02_control_ifg_run%d_timepoints.txt',basedir,sub_current));
        pmtg_tmp=load(sprintf('%/mps/stand/roiTP/sub%02_control_pmtg_run%d_timepoints.txt',basedir,sub_current));
        rpmtg_tmp=load(sprintf('%/mps/stand/roiTP/sub%02_pmtg_seman_run%d_timepoints.txt',basedir,sub_current));
        ifg=[ifg;ifg_tmp];
        pmtg=[pmtg;pmtg_tmp];
        rpmtg=[rpmtg;rpmtg_tmp];
    end
    scact=[ifg pmtg rpmtg];
    
    
    filename_ps = sprintf('%s/mps/ptRSA/tgroup/sub%02d_allrsa.nii.gz',basedir,sub_current);
    allrsa=niftiread(filename_ps);
    
    filename_act = sprintf('%s/mps/ptRSA/tgroup/sub%02d_allact.nii.gz',basedir,sub_current);
    allact=niftiread(filename_act);
    
    
    allrsa=shiftdim(allrsa,3);
    allact=shiftdim(allact,3);


    allsta=zeros(xlength,ylength,zlength,9);
    
    
    
    for k=radius+1:step:xlength-radius
        for j=radius+1:step:ylength-radius
            for i=radius+1:step:zlength-radius
                
                cubic_p       = allrsa(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_p = cubic_p(:,:);
                mps=mean(cubicvector_p');
                
                cubic_a       = allact(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_a = cubic_a(:,:);
                mact=mean(cubicvector_a');
                
                                
                
                if min(std(cubicvector_p,0,2))<epsilon || min(std(cubicvector_a,0,2))<epsilon
                    
                    allsta(k,j,i,:)=10;
                    
                else
                    
                    
                    for r=1:3
                        allsta(k,j,i,(r-1)*3+1)=partialcorr(mps(alldata.KeyPress==1),scact(alldata.KeyPress==1,r),mact(alldata.KeyPress==1));
                        allsta(k,j,i,(r-1)*3+2)=partialcorr(mps(alldata.KeyPress==2),scact(alldata.KeyPress==2,r),mact(alldata.KeyPress==2));
                        allsta(k,j,i,(r-1)*+3)=partialcorr(mps,scact(:,r),mact);
                    end
                    
                    
                    
                    
                    
                    
                end % end
            end % end
        end % end
    end
    
    label={'control_ifg_y','control_ifg_n','control_ifg_a','pmtg_y','pmtg_n','pmtg_a','control_pmtg_y','control_pmtg_n','control_pmtg_a'};
    for l=1:length(label)
        filename=sprintf('%s/mps/controlPS/tgroup/sub%02d_%s.nii',basedir,sub_current,label{l});
        
        data=squeeze(allrsa(:,:,:,l));
        niftiwrite(data,filename)
        system(sprintf('gzip -f %s',filename));
    end
end
end
