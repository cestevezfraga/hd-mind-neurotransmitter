
%Script for statistics 12/03/14   
    
load 'group_results'
load 'age'
load 'sex'
age = transpose(age);
sex = transpose(sex);

%Permutation test for degree
design=[ones(1,17) 2*ones(1,18)];
N=length(design);
T_deg=cat(3,Tcont_deg,Tpat_deg);
T_deg=reshape(T_deg,80,35); %all subj degree
data=T_deg;
for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[age' ones(length(age),1)]);
    data(n,:)=R; % we keep only the residual
end

for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[sex' ones(length(sex),1)]);
    data(n,:)=R; % we keep only the residual
end
niter=10000;
	stats.deg = bramila_ttest2_np(data,design,niter);
	p_fdr_deg=mafdr(stats.deg.pvals(:,1));
	p_fdr_deg(p_fdr_deg>0.05)=0;

%Permutation test for strength
T_stren=cat(3,Tcont_stren,Tpat_stren);
T_stren=reshape(T_stren,80,35);
data=T_stren;
for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[age' ones(length(age),1)]);
    data(n,:)=R; % we keep only the residual
end

for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[sex' ones(length(sex),1)]);
    data(n,:)=R; % we keep only the residual
end
niter=10000;
	stats.stren = bramila_ttest2_np(data,design,niter);
	p_fdr_stren=mafdr(stats.stren.pvals(:,1));
	p_fdr_stren(p_fdr_stren>0.05)=0;
    
% Permutation test for local efficiency    
T_Eloc=cat(3,Tcont_Eloc,Tpat_Eloc);
T_Eloc=reshape(T_Eloc,80,35);
data=T_Eloc;
for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[age' ones(length(age),1)]);
    data(n,:)=R; % we keep only the residual
end

for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[sex' ones(length(sex),1)]);
    data(n,:)=R; % we keep only the residual
end
niter=10000;
	stats.Eloc = bramila_ttest2_np(data,design,niter);
	p_fdr_Eloc=mafdr(stats.Eloc.pvals(:,1));
	p_fdr_Eloc(p_fdr_Eloc>0.05)=0;
	 
%Permutation test for coefficient clustering
T_CC=cat(3,Tcont_CC,Tpat_CC);
T_CC=reshape(T_CC,80,35);
data=T_CC;
for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[age' ones(length(age),1)]);
    data(n,:)=R; % we keep only the residual
end

for n=1:size(data,1)    % for each data point
    y=data(n,:);
    [B,BINT,R]=regress(y',[sex' ones(length(sex),1)]);
    data(n,:)=R; % we keep only the residual
end
niter=10000;
	stats.CC_P = bramila_ttest2_np(data,design,niter);
	p_fdr_CC=mafdr(stats.CC_P.pvals(:,1));
	p_fdr_CC(p_fdr_CC>0.05)=0;
	
Stats_results=cat(2,p_fdr_deg,p_fdr_stren,p_fdr_Eloc,p_fdr_CC);
