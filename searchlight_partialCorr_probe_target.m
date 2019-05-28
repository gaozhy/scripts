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
    
    load(sprintf('%/roiAct/sem/sub%02_controlAct.mat',basedir,sub_current));
    
    
    
    
    filename_ps = sprintf('%s/mps/ptRSA/tgroup/sub%02d_allrsa.nii.gz',basedir,sub_current);
    allrsa=niftiread(filename_ps);
    
    filename_std = sprintf('%s/mps/ptRSA/tgroup/sub%02d_allact.nii.gz',basedir,sub_current);
    allact=niftiread(filename_std);
    
    filename_std = sprintf('%s/mps/ptRSA/tgroup/sub%02d_allstd.nii.gz',basedir,sub_current);
    allstd=niftiread(filename_std);
    
    
    allrsa=shiftdim(allrsa,3);
    allact=shiftdim(allact,3);
    allstd=shiftdim(allstd,3);

    allsta=zeros(xlength,ylength,zlength,7);
    
    
    
    for k=radius+1:step:xlength-radius
        for j=radius+1:step:ylength-radius
            for i=radius+1:step:zlength-radius
                
                cubic_p       = allrsa(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_p = cubic_p(:,:);
                cubic1=mean(cubicvector_p);
                tmp1=cubicvector_p';
                
                cubic_t       = allact(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_t = cubic_t(:,:);
                cubic2=mean(cubicvector_t);
                tmp2=cubicvector_t';
                
                cubic_d       = allstd(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_d = cubic_d(:,:);
                cubic3=mean(cubicvector_d);
                tmp3=cubicvector_d';
                                
                
                if min(std(cubicvector_p,0,2))<epsilon || min(std(cubicvector_t,0,2))<epsilon 

                    
                    allsta(k,j,i,:)=10;
                    
                else
                    cc=diag(corr(cubicvector_p',cubicvector_t'));
                    z_cc=0.5*(log(1+cc)-log(1-cc));
                    
                    
    
                    act=mean(tmp1+tmp2)/2;act=act';%std=nanmean(nanstd(tmp1)+nanstd(tmp2))/2;
                    allstd(k,j,i,:)=nanmean(nanstd(tmp1)+nanstd(tmp2))/2;
                    
                    corr_y=partialcorr(z_cc(alldata.KeyPress==1),allwb(alldata.KeyPress==1),act(alldata.KeyPress==1));
                    corr_n=partialcorr(z_cc(alldata.KeyPress==2),allwb(alldata.KeyPress==2),act(alldata.KeyPress==2));
                    corr_a=partialcorr(z_cc,allwb(:,3),act);
                    
                    zcorr_y=0.5*(log(1+corr_y)-log(1-corr_y));
                    zcorr_n=0.5*(log(1+corr_n)-log(1-corr_n));
                    zcorr_a=0.5*(log(1+corr_a)-log(1-corr_a));
                    
                    
                    
                    allsta(k,j,i,1)=mean(z_cc(alldata.KeyPress==1));
                    allsta(k,j,i,2)=mean(z_cc(alldata.KeyPress==2));
                    allsta(k,j,i,3)=mean(z_cc(alldata.KeyPress==1))-mean(z_cc(alldata.KeyPress==2));
                    allsta(k,j,i,4)=zcorr_a;
                    allsta(k,j,i,5)=zcorr_y;
                    allsta(k,j,i,6)=zcorr_n;
                    allsta(k,j,i,7)=zcorr_y-zcorr_n;
                    
                    
                    
                end % end
            end % end
        end % end
    end
    
    label={'ps_y','ps_n','ps_ymn','corr_a','corr_y','corr_n','corr_ymn'};
    for l=1:length(label)
        filename=sprintf('%s/mps/ptRSA/ptgroup/sub%02d_%s.nii',basedir,sub_current,label{l});
        
        data=squeeze(allrsa(:,:,:,l));
        niftiwrite(data,filename)
        system(sprintf('gzip -f %s',filename));
    end
    filename1=sprintf('%s/mps/ptRSA/ptgroup/sub%02d_allrsa.nii',basedir,sub_current);
    niftiwrite(allrsa,filename1);
    system(sprintf('gzip -f %s',filename1));
    
    filename2=sprintf('%s/mps/ptRSA/ptgroup/sub%02d_allact.nii',basedir,sub_current);
    niftiwrite(allact,filename2);
    system(sprintf('gzip -f %s',filename2));
    
    filename3=sprintf('%s/mps/ptRSA/ptgroup/sub%02d_allstd.nii',basedir,sub_current);
    niftiwrite(allstd,filename3);
    system(sprintf('gzip -f %s',filename3));
    
    time=toc/3600;
end
end
