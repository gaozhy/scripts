function searchlight_Beh_Brain_CP_probe_target
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


load(sprintf('%s/beh_results/sem/allw2b.mat',basedir)); %load word-pairs' semantic association strength and level information, col-1:trial_id, col-2:w2b, col-3~8: num of levels 3~8; col-9: runid, col-10: pres_id when exp running
load(sprintf('%s/beh_results/sem/w2b_beh_struc.mat',basedir)); % load behavioral level similarity structure for probe and target separately.
load(sprintf('%s/beh_results/sem/biased_beh_struc.mat',basedir)); % load biased beh-similarity structure


for iSub=1:length(subs)
    sub_current=subs(iSub);
    fprintf('Subject %d is processing.\n',sub_current);
    
    load(sprintf('%s/beh_results/sem/sub%02d_alldata.mat',basedir,sub_current));
    
    weak_no=alldata.KeyPress==2&allwb(:,4)==1;
    
    
    
    
    filename_probe = sprintf('%s/mps/stand/sem/sub%02d_probe_beta_allruns.nii.gz',basedir,sub_current);
    all_pdata=niftiread(filename_probe);
    filename_target = sprintf('%s/mps/stand/sem/sub%02d_target_beta_allruns.nii.gz',basedir,sub_current);
    all_tdata=niftiread(filename_target);
    
    
    all_pdata=shiftdim(all_pdata,3);% do matrix shift, so the trial list become the first dimension;
    all_tdata=shiftdim(all_tdata,3);
    

    allsta=zeros(xlength,ylength,zlength,7);
    
    
    
    for k=radius+1:step:xlength-radius
        for j=radius+1:step:ylength-radius
            for i=radius+1:step:zlength-radius
                
                cubic_p       = all_pdata(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_p = cubic_p(:,:);
                cubic1=mean(cubicvector_p);
                tmp1=cubicvector_p';
                
                cubic_t       = all_tdata(:,k-radius:k+radius,j-radius:j+radius,i-radius:i+radius); % define small cubic
                cubicvector_t = cubic_t(:,:);
                cubic2=mean(cubicvector_t);
                tmp2=cubicvector_t';
                
                if max(abs(cubic1-mean(cubic1)))<epsilon || max(cubic2-mean(cubic2))<epsilon

                    allsta(k,j,i,:)=10000;
                else
                    cc_t=corr(tmp2,tmp2);
                    cc_p=corr(tmp1,tmp1);
                    
                    z_ccp=0.5*(log(1+cc_p)-log(1-cc_p));
                    z_cct=0.5*(log(1+cc_t)-log(1-cc_t));
                    
                    weak_nop=[];
                    weak_not=[];
                    weak_no_bp=[];
                    weak_no_bt=[];
                    weak_no_bb=[];
                    
                    medst_l1p=[];
                    medst_l1t=[];                 
                    medst_l1_bp=[];
                    medst_l1_bt=[];
                    medst_l1_bb=[];
                    
                    medst_l2p=[];
                    medst_l2t=[];                 
                    medst_l2_bp=[];
                    medst_l2_bt=[];
                    medst_l2_bb=[];
                    
                    medst_l3p=[];
                    medst_l3t=[];                 
                    medst_l3_bp=[];
                    medst_l3_bt=[];
                    medst_l3_bb=[];                    
                    
                    
                    for r=1:3
                        
                        
                        weak_no_c=alldata.KeyPress==2&allwb(:,4)==1&alldata.runid==r;
                        weak_no_r=alldata.KeyPress==2&allwb(:,4)==1&alldata.runid>r;
                        medst_l1_c=alldata.KeyPress==1&allwb(:,4)==2&alldata.runid==r;
                        medst_l1_r=alldata.KeyPress==1&allwb(:,4)==2&alldata.runid>r;
                        medst_l2_c=alldata.KeyPress==1&allwb(:,4)==3&alldata.runid==r;
                        medst_l2_r=alldata.KeyPress==1&allwb(:,4)==3&alldata.runid>r;
                        medst_l3_c=alldata.KeyPress==1&allwb(:,4)==4&alldata.runid==r;
                        medst_l3_r=alldata.KeyPress==1&allwb(:,4)==4&alldata.runid>r;
                        
                        weak_nop_tmp=z_ccp(weak_no_c,weak_no_r);
                        medst_l1p_tmp=z_ccp(medst_l1_c,medst_l1_r);
                        medst_l2p_tmp=z_ccp(medst_l2_c,medst_l2_r);
                        medst_l3p_tmp=z_ccp(medst_l3_c,medst_l3_r);                        
                        
                        weak_not_tmp=z_cct(weak_no_c,weak_no_r);
                        medst_l1t_tmp=z_cct(medst_l1_c,medst_l1_r);
                        medst_l2t_tmp=z_cct(medst_l2_c,medst_l2_r);
                        medst_l3t_tmp=z_cct(medst_l3_c,medst_l3_r);

                        
                        
                        weak_no_bp_tmp=w2b_struc_p(weak_no_c,weak_no_r);
                        medst_l1_bp_tmp=w2b_struc_p(medst_l1_c,medst_l1_r);
                        medst_l2_bp_tmp=w2b_struc_p(medst_l2_c,medst_l2_r);
                        medst_l3_bp_tmp=w2b_struc_p(medst_l3_c,medst_l3_r);                        
                        
                        weak_no_bt_tmp=w2b_struc_t(weak_no_c,weak_no_r);
                        medst_l1_bt_tmp=w2b_struc_t(medst_l1_c,medst_l1_r);
                        medst_l2_bt_tmp=w2b_struc_t(medst_l2_c,medst_l2_r);
                        medst_l3_bt_tmp=w2b_struc_t(medst_l3_c,medst_l3_r);                       
                        
                        weak_no_bb_tmp=w2b_struc_b(weak_no_c,weak_no_r);
                        medst_l1_bb_tmp=w2b_struc_b(medst_l1_c,medst_l1_r);
                        medst_l2_bb_tmp=w2b_struc_b(medst_l2_c,medst_l2_r);
                        medst_l3_bb_tmp=w2b_struc_b(medst_l3_c,medst_l3_r); 
                        
                        weak_nop=[weak_nop;weak_nop_tmp];
                        weak_not=[weak_not;weak_not_tmp];
                        
                        weak_no_bp=[weak_no_bp;weak_no_bp_tmp];
                        weak_no_bt=[weak_no_bt;weak_no_bt_tmp];
                        weak_no_bb=[weak_no_bb;weak_no_bb_tmp];
                        
                        
                        medst_l1p=[medst_l1p;medst_l1p_tmp];
                        medst_l1t=[medst_l1t;medst_l1t_tmp];
                        
                        medst_l1_bp=[medst_l1_bp;medst_l1_bp_tmp];
                        medst_l1_bt=[medst_l1_bt;medst_l1_bt_tmp];                    
                        medst_l1_bb=[medst_l1_bb;medst_l1_bb_tmp];
                        
                        
                        
                        medst_l2p=[medst_l2p;medst_l2p_tmp];
                        medst_l2t=[medst_l2t;medst_l2t_tmp];
                        
                        medst_l2_bp=[medst_l2_bp;medst_l2_bp_tmp];
                        medst_l2_bt=[medst_l2_bt;medst_l2_bt_tmp];                    
                        medst_l2_bb=[medst_l2_bb;medst_l2_bb_tmp];
                        
                        
                        
                        medst_l3p=[medst_l3p;medst_l3p_tmp];
                        medst_l3t=[medst_l3t;medst_l3t_tmp];
                        
                        medst_l3_bp=[medst_l3_bp;medst_l3_bp_tmp];
                        medst_l3_bt=[medst_l3_bt;medst_l3_bt_tmp];                    
                        medst_l3_bb=[medst_l3_bb;medst_l3_bb_tmp];                        
                        
                        
                    end
                    


                    
                    
                    
                    weak_no_corr_tp=corr(weak_not,weak_no_bp);
                    weak_no_corr_tt=corr(weak_not,weak_no_bt);
                    weak_no_corr_tb=corr(weak_not,weak_no_bb);
                    
                    medst_l1_corr_tp=corr(medst_l1t,medst_l1_bp);
                    medst_l1_corr_tt=corr(medst_l1t,medst_l1_bt);
                    medst_l1_corr_tb=corr(medst_l1t,medst_l1_bb);
                    
                    medst_l2_corr_tp=corr(medst_l2t,medst_l2_bp);
                    medst_l2_corr_tt=corr(medst_l2t,medst_l2_bt);
                    medst_l2_corr_tb=corr(medst_l2t,medst_l2_bb);
                    
                    medst_l3_corr_tp=corr(medst_l3t,medst_l3_bp);
                    medst_l3_corr_tt=corr(medst_l3t,medst_l3_bt);
                    medst_l3_corr_tb=corr(medst_l3t,medst_l3_bb);
                    
                    
                    weak_no_corr_pp=corr(weak_nop,weak_no_bp);
                    weak_no_corr_pt=corr(weak_nop,weak_no_bt);
                    weak_no_corr_pb=corr(weak_nop,weak_no_bb);
                    
                    
                    medst_l1_corr_pp=corr(medst_l1p,medst_l1_bp);
                    medst_l1_corr_pt=corr(medst_l1p,medst_l1_bt);
                    medst_l1_corr_pb=corr(medst_l1p,medst_l1_bb);
                    
                    medst_l2_corr_pp=corr(medst_l2p,medst_l2_bp);
                    medst_l2_corr_pt=corr(medst_l2p,medst_l2_bt);
                    medst_l2_corr_pb=corr(medst_l2p,medst_l2_bb);
                    
                    medst_l3_corr_pp=corr(medst_l3p,medst_l3_bp);
                    medst_l3_corr_pt=corr(medst_l3p,medst_l3_bt);
                    medst_l3_corr_pb=corr(medst_l3p,medst_l3_bb);
                    
                    
                    
                    
     
                    
                    allsta(k,j,i,1)=weak_no_corr_tp;
                    allsta(k,j,i,2)=weak_no_corr_tt;
                    allsta(k,j,i,3)=weak_no_corr_tb;
                    allsta(k,j,i,4)=medst_l1_corr_tp;
                    allsta(k,j,i,5)=medst_l1_corr_tt;
                    allsta(k,j,i,6)=medst_l1_corr_tb;
                    allsta(k,j,i,7)=medst_l2_corr_tp;
                    allsta(k,j,i,8)=medst_l2_corr_tt;
                    allsta(k,j,i,9)=medst_l2_corr_tb;
                    allsta(k,j,i,10)=medst_l3_corr_tp;
                    allsta(k,j,i,11)=medst_l3_corr_tt;
                    allsta(k,j,i,12)=medst_l3_corr_tb;
                    
                    allsta(k,j,i,13)=weak_no_corr_pp;
                    allsta(k,j,i,14)=weak_no_corr_pt;
                    allsta(k,j,i,15)=weak_no_corr_pb;
                    allsta(k,j,i,16)=medst_l1_corr_pp;
                    allsta(k,j,i,17)=medst_l1_corr_pt;
                    allsta(k,j,i,18)=medst_l1_corr_pb;
                    allsta(k,j,i,19)=medst_l2_corr_pp;
                    allsta(k,j,i,20)=medst_l2_corr_pt;
                    allsta(k,j,i,21)=medst_l2_corr_pb;
                    allsta(k,j,i,22)=medst_l3_corr_pp;
                    allsta(k,j,i,23)=medst_l3_corr_pt;
                    allsta(k,j,i,24)=medst_l3_corr_pb;                    

                    
                    
                end % end
            end % end
        end % end
    end
    
    label={'weak_no_corr_tp','weak_no_corr_tt','weak_no_corr_tb','medst_l1_corr_tp','medst_l1_corr_tt','medst_l1_corr_tb',...
        'medst_l2_corr_tp','medst_l2_corr_tt','medst_l2_corr_tb','medst_l3_corr_tp','medst_l3_corr_tt','medst_l3_corr_tb',...
        'weak_no_corr_pp','weak_no_corr_pt','weak_no_corr_pb','medst_l1_corr_pp','medst_l1_corr_pt','medst_l1_corr_pb',...
        'medst_l2_corr_pp','medst_l2_corr_pt','medst_l2_corr_pb','medst_l3_corr_pp','medst_l3_corr_pt','medst_l3_corr_pb'};
    
    for l=1:length(label)
        filename=sprintf('%s/mps/beh_brain_constrain/bgroup/sub%02d_%s.nii',basedir,sub_current,label{l});
        
        data=squeeze(allsta(:,:,:,l));
        niftiwrite(data,filename)
        system(sprintf('gzip -f %s',filename));
    end

    
    time=toc/3600;
end
end