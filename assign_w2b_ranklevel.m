clear;clc;

basedir='C:\Users\zhiya\Documents\MATLAB\beh_results'
beh_dir=fullfile(basedir,'sem');

cd(beh_dir);

load('sub01_alldata.mat');

w2b=alldata.word2vec;
trial_id=alldata.Trial_id;
pres_id=[1:192]';
allwb=zeros(192,10);
allwb(:,1)=trial_id;
allwb(:,2)=w2b;
allwb(:,10)=pres_id;

allwb_sorted=sortrows(allwb,2);

num_trial_level=64;
num_levels=192/num_trial_level;

for i=1:num_levels
    allwb_sorted((i-1)*num_trial_level+1:i*num_trial_level,3)=i;
end

num_trial_level=48;

num_levels=192/num_trial_level;

for i=1:num_levels
    allwb_sorted((i-1)*num_trial_level+1:i*num_trial_level,4)=i;
end

num_trial_level=38;

num_levels=192/num_trial_level;

for i=1:num_levels
    a=mod(192,num_trial_level);
    if a~=0
    allwb_sorted((i-1)*num_trial_level+3:i*num_trial_level+2,5)=i;
    end
end

num_trial_level=32;

num_levels=192/num_trial_level;

for i=1:num_levels
    allwb_sorted((i-1)*num_trial_level+1:i*num_trial_level,6)=i;
end

num_trial_level=27;

num_levels=192/num_trial_level;

for i=1:num_levels
    a=mod(192,num_trial_level);
    if a~=0
    allwb_sorted((i-1)*num_trial_level+3:i*num_trial_level+2,7)=i;
    end
end
num_trial_level=24;

num_levels=192/num_trial_level;

for i=1:num_levels
    allwb_sorted((i-1)*num_trial_level+1:i*num_trial_level,8)=i;
end

allwb=sortrows(allwb_sorted,10);
for i=1:4
    allwb((i-1)*48+1:i*48,9)=i;
end