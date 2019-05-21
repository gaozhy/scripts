clear;clc;

basedir='C:\Users\zhiya\Documents\MATLAB\beh_results'
beh_dir=fullfile(basedir,'wm');

cd(beh_dir);

subs=[1:8 10:20];

levels=[3,4,5,6,7];

allrt_y=zeros(length(subs),5);allratio_y=zeros(length(subs),5);
allrt_n=zeros(length(subs),5);allratio_n=zeros(length(subs),5);

allrt_hit=zeros(length(subs),5);
allrt_cr=zeros(length(subs),5);

allratio_hit=zeros(length(subs),5);
allratio_cr=zeros(length(subs),5);


for l=3:7
    

    
  

    

       
        
        for sub=1:length(subs)
            subid=subs(sub);
            filename=sprintf('sub%02d_alldata.mat',subid);
            load(filename);
            
            
            
            
            rt_y_tmp=mean(alldata.RT(alldata.Level==l&alldata.KeyPress==alldata.Correct_Key));
            rt_n_tmp=mean(alldata.RT(alldata.Level==l&alldata.KeyPress~=alldata.Correct_Key));
            
            ratio_y_tmp=sum(alldata.Level==l&alldata.KeyPress==alldata.Correct_Key)/length(find(alldata.Level==l));
            ratio_n_tmp=sum(alldata.Level==l&alldata.KeyPress~=alldata.Correct_Key)/length(find(alldata.Level==l));
            
            allrt_hit(sub,l-2)=mean(alldata.RT(alldata.Level==l&alldata.KeyPress==1&alldata.Correct_Key==1));
            allrt_cr(sub,l-2)=mean(alldata.RT(alldata.Level==l&alldata.KeyPress==2&alldata.Correct_Key==2));
            
            allrt_y(sub,l-2)=rt_y_tmp;
            allrt_n(sub,l-2)=rt_n_tmp;
            
            allratio_y(sub,l-2)=ratio_y_tmp;
            allratio_n(sub,l-2)=ratio_n_tmp;
            
            allratio_cr(sub,l-2)=sum(alldata.Level==l&alldata.KeyPress==2&alldata.Correct_Key==2)/length(find(alldata.Level==l&alldata.Correct_Key==2));
            allratio_hit(sub,l-2)=sum(alldata.Level==l&alldata.KeyPress==1&alldata.Correct_Key==1)/length(find(alldata.Level==l&alldata.Correct_Key==1));
            
            
            
        end
        
        
   

end
    filename_out=sprintf('beh_results_wm_levels.mat');
    eval(sprintf('save %s allrt_y allrt_n allratio_y allratio_n allrt_hit allrt_cr allratio_hit allratio_cr',filename_out));