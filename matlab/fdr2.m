function [varargout]=fdr2(pval,varargin)
%% FDR false discovery rate correction for multiple comparison.
%
% USAGE: [FDR_threshold,survived,uncorreced_p_for_survived,...
%          adjusted_p_for_survived, all_adjusted_p]=fdr(p_values, [q,cv]);
%
%         Here, p-value could be from any of multiple comparison test. 
%         However p-value should be either a vector or [edges,measures].
%
% False discovery rate correction , Written by Cheol E Han on Feb 9, 2011
%                           Last Modified by Cheol E Han, on Feb 26, 2012
% Related article: 
%   Benjamini Y and Hochberg Y (1995): Controlling the False Discovery Rate:
%      A Practical and Powerful Approach to Multiple Testing. J R Statist 
%      Soc B. 57(1):289-300.
%   Yekutieli D and Benjamini Y (1999): Resampling-based false discovery 
%      rate controlling multiple test procedures for correlated test stati-
%      -stics. J Stat Plan Infer. 82(1-2):171-96.
%   Benjamini Y and Yekutieli D (2001): The control of the false discovery 
%      rate in multiple testing under dependency. Ann Statist. 29(4):1165-88.
%   Genovese CR, Lazar NA, Nichols T (2002): Thresholding of statistical 
%      maps in functional neuroimaging using the false discovery rate. 
%      Neuroimage 15:870 ? 878.
% 
% See also ttest2, perm_test
% 
%


%% Function body
if nargout==0;	error('no output arguments');	end
if numel(pval)==length(pval) && size(pval,1)<size(pval,2);     pval=pval';  end
n=size(pval,1);
no_measures=size(pval,2);

% at least two arguments are needed
qval=[]; cv=[];
switch nargin
    case 1                  % cv for independency
        qval=0.05;
        cv=sum(1./(1:n));    
    case 2                  % cv for dependency (default)
        qval=varargin{1};      % following Benjamini Y and Yekutieli D (2001)
        cv=sum(1./(1:n));   % and Genovese CR, Lazar NA, Nichols T (2002)
    case 3
        qval=varargin{1};
        if varargin{2}==1;  cv=1;     end   
end
if numel(qval)==0;      qval=0.05;              end
if numel(cv)==0;        cv=sum(1./(1:n));       end
if qval>1 || qval<0;          
    qval=0.05;   
    fprintf('q-val should be between 0 and 1. set to a default value=0.05. (0.2 is in genreal maximum.)\n');  
end

thr=zeros(no_measures,1);
survived_set=cell(no_measures,1);
p_of_survived_set=cell(no_measures,1);
ap_of_survived_set=cell(no_measures,1);
p_adjusted=zeros(n,no_measures);

% fdr related parameter 
tol=(1:n)'/n*qval/cv;
for m=1:no_measures
    [thr(m),survived,survived_p,p_adjusted(:,m)]=fdr_p(pval(:,m),tol,cv,n);
    survived_set{m}=survived;
    p_of_survived_set{m}=survived_p;
    ap_of_survived_set{m}=p_adjusted(survived,m);
end

%% FDR_threhold: maximum p in uncorrected p-value 
% This result is the same as outputs of Nichols' fdr outine (FDR_nichols) 
% This result is the same as FDR_threhold of Anderson M. Winkler's routine (FDR_aw)
if nargout>=1;   varargout(1)={thr};   end

%% survived edges, ordering with smallest uncorreted p-values.
if nargout>=2;   
    if no_measures==1;   varargout(2)=survived_set(1);   
    else varargout(2)={survived_set};   end
end

%% uncorrected p-values for survived edges
if nargout>=3;
    if no_measures==1;   varargout(3)=p_of_survived_set(1);
    else varargout(3)={p_of_survived_set};                   end
end

%% adjusted p-value for only survived edges
if nargout>=4;
    if no_measures==1;   varargout(4)=ap_of_survived_set(1);
    else varargout(4)={ap_of_survived_set};                  end
end

%% adjusted p-value for all edges
if nargout>=5;  varargout(5)={p_adjusted};               end





%----------------------private function-----------------------------------
function [thr,survived,p_of_survived,p_adjusted]=fdr_p(p,tol,cv,n)

[sorted,sort_idx]=sort(p);
max_p_idx=find(sorted<=tol,1,'last');
thr=sorted(max_p_idx);
survived=sort_idx(1:max_p_idx);
p_of_survived=sorted(1:max_p_idx);
if numel(max_p_idx)==0;     thr=0;      end

%% FDR adjusted p-value
% following Anderson M. Winkler's webpage (see below)
% FDR corrected p-value is indeed incorrect one.
% FDR adjusted p-value is more sensible and correct. 
% (http://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/)
%
% This follows Yekutieli & Benjamini (1999), equation #3 and
% code in http://www0.cs.ucl.ac.uk/staff/gridgway/stats/fdr_query.html
Qs=sorted*n./(1:n)'*cv;
p_adjusted    = ones(size(p));
for i=1:n
    idx=find(sorted>=p(i),1,'first');
    if ~isempty(idx);    p_adjusted(i)=min(Qs(idx:n));     end
end