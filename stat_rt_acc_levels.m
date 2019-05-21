clear;clc;

basedir='C:\Users\zhiya\Documents\MATLAB\beh_results'
beh_dir=fullfile(basedir,'sem');

cd(beh_dir);

subs=1:20;
load('allw2b.mat');
levels=[4,5,6,7,8];


for l=1:length(levels)
    num_level=levels(l);
    allrt_y=zeros(length(subs),num_level);allratio_y=zeros(length(subs),num_level);
    allrt_n=zeros(length(subs),num_level);allratio_n=zeros(length(subs),num_level);
    
  

    
    for i=1:num_level
       
        
        for sub=1:length(subs)
            subid=subs(sub);
            filename=sprintf('sub%02d_alldata.mat',subid);
            load(filename);
            
            
            
            
            rt_y_tmp=mean(alldata.RT(allwb(:,num_level)==i&alldata.KeyPress==1));
            rt_n_tmp=mean(alldata.RT(allwb(:,num_level)==i&alldata.KeyPress==2));
            
            ratio_y_tmp=sum(allwb(:,num_level)==i&alldata.KeyPress==1)/length(find(allwb(:,num_level)==i));
            ratio_n_tmp=sum(allwb(:,num_level)==i&alldata.KeyPress==2)/length(find(allwb(:,num_level)==i));
            
            
            allrt_y(sub,i)=rt_y_tmp;
            allrt_n(sub,i)=rt_n_tmp;
            
            allratio_y(sub,i)=ratio_y_tmp;
            allratio_n(sub,i)=ratio_n_tmp;
            
            
            
            
        end
        
        
    end
    filename_out=sprintf('beh_results_%dlevel.mat',num_level);
    eval(sprintf('save %s allrt_y allrt_n allratio_y allratio_n',filename_out));
end
