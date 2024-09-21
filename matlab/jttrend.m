function jttrend(varargin)
% JTTREND: Perform the Jonckheere-Terpstra test on trend.
% There are situations in which treatments are ordered in some
% way, for example the increasing dosages of a drug. In these
% cases a test with the more specific alternative hypothesis that
% the population medians are ordered in a particular direction
% may be required. For example, the alternative hypothesis
% could be as follows: population median1 <= population
% median2 <= population median3. This is a one-tail test, and
% reversing the inequalities gives an analagous test in the
% opposite tail. Here, the Jonckheere–Terpstra test can be
% used.
% Bewick V., Cheek L., Ball J. Statistics review 10: further nonparametric
% methods. Critical Care 2004, 8: 196-199
%
% Assumptions:
% - Data must be at least ordinal
% - Groups must be selected in a meaningful order i.e. ordered
% If you do not choose to enter your own group scores then scores are
% allocated uniformly (1 ... n) in order of selection of the n groups.
%
% Syntax: 	jttrend(x,alpha,score)
%      
%     Inputs:
%           X - Nx2 data matrix 
%           ALPHA - significance level (default 0.05) 
%           SCORE - order of selection of the groups
%     Outputs:
%           - Jonckheere-Terpstra statistics and p-value
%
%   Example:
% Mice were inoculated with cell lines, CMT 64 to 181, which had been
% selected for their increasing metastatic potential. The number of lung
% metastases found in each mouse after inoculation are quoted below:
%
%                                 Sample
%                   ---------------------------------
%                      64   167  170  175  181
%                   ---------------------------------
%                      0    0    2    0    2
%                      0    0    3    3    4
%                      1    5    6    5    6
%                      1    7    9    6    6
%                      2    8    10   10   6
%                      2    11   11   19   7
%                      4    13   11   56   18
%                      9    23   12   100  39    
%                           25   21   132  60
%                           97
%                   ---------------------------------
%
%       Data matrix must be:
%    d=[0 0 1 1 2 2 4 9 0 0 5 7 8 11 13 23 25 97 2 3 6 9 10 11 11 12 21 ...
%       0 3 5 6 10 19 56 100 132 2 4 6 6 6 7 18 39 60];
%    g=[ones(1,8) 2.*ones(1,10) 3.*ones(1,9) 4.*ones(1,9) 5.*ones(1,9)];
%    x=[d' g'];
%
%           Calling on Matlab the function: jttrend(x)
% (in this case, the groups are automated scored from 1 to 5)
%
%           Answer is:
%
% JONCKHEERE-TERPSTRSA TEST FOR NON PARAMETRIC TREND ANALYSIS
%  
% Comparison			Nx		Ny			Uxy
% 1-2					8		10			63.0
% 1-3					8		9			65.5
% 1-4					8		9			61.0
% 1-5					8		9			63.5
% 2-3					10		9			41.0
% 2-4					10		9			49.5
% 2-5					10		9			41.5
% 3-4					9		9			45.5
% 3-5					9		9			39.0
% 4-5					9		9			35.5
% 										    ----- 
% 								       U = 505.0
%  
% JT			1-sided p			H0
% 2.0116		 0.0221			 rejected.
% We have shown a statistically significant trend for increasing number of
% metastases across these malignant cell lines in this order.
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Jonckheere-Terpstra test: A nonparametric Test for Trend
% http://www.mathworks.com/matlabcentral/fileexchange/22159

%Input Error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu)
    error('Warning: almost data matrix is required')
elseif nu>3
    error('Warning: Max three input data are required')
end
default.values = {[],0.05,[]};
default.values(1:nu) = args;
[x alpha score] = deal(default.values{:});

if isempty(x)
    error('Warning: X matrix is empty...')
end
if isvector(x)
    error('Warning: x must be a matrix, not a vector.');
end
if ~all(isfinite(x(:))) || ~all(isnumeric(x(:)))
    error('Warning: all X values must be numeric and finite')
end
%check if x is a Nx2 matrix
if size(x,2) ~= 2
    error('Warning: JTTREND requires a Nx2 input matrix')
end
%check if x(:,2) are all whole elements
if ~all(x(:,2) == round(x(:,2)))
    error('Warning: all elements of column 2 of input matrix must be whole numbers')
end
%check if there are 3 groups at least
k=max(x(:,2)); %number of groups
if k<3
    error('Warning: Three groups at least are required')
end
if nu>=2
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
if nu==3 && ~isempty(score) %check score
    if ~isvector(score) || ~all(isfinite(score)) || ~all(isnumeric(score)) || ~all(score == round(score))
        error('Warning: SCORE is a vector and all its values must be numeric, finite and integer')
    end
else
    score=1:1:k;
end
clear args default nu

ni=crosstab(x(:,2)); %elements for each group
N=sum(ni); %total observation
%Build-up the matrix of observation
X=repmat(NaN,max(ni),k);
for I=1:k
    X(1:ni(score(I)),I)=x(x(:,2)==score(I),1);
end
%vector and variable preallocation
comp=0.5*k*(k-1); %number of comparison
U=zeros(1,comp); G=1;

disp('JONCKHEERE-TERPSTRSA TEST FOR NON PARAMETRIC TREND ANALYSIS')
disp(' ')
fprintf('Comparison\t\t\tNx\t\tNy\t\t\tUxy\n')
for I=1:k-1
    Uxy=zeros(1,ni(score(I)));
    for J=I+1:k
        for F=1:ni(score(I)) %for each element of the i-esim group...
            %...check how many elements of j-esim group is major or equal
            %then (in this case multiple them for 1/2)
            Uxy(F)=length(X(X(:,J)>X(F,I)))+0.5*length(X(X(:,J)==X(F,I)));
        end
        U(G)=sum(Uxy); 
        fprintf('%d-%d\t\t\t\t\t%d\t\t%d\t\t\t%0.1f\n',I,J,ni(score(I)),ni(score(J)),U(G))
        G=G+1;
    end
end
Ut=sum(U); N2=N^2; ni2=ni.^2;
fprintf('\t\t\t\t\t\t\t\t\t\t----- \n')
fprintf('\t\t\t\t\t\t\t\t    U = %0.1f\n',Ut)
disp(' ')
%Compute the Jonckheere-Terpstra statistics
num=Ut-(N2-sum(ni2))/4;
denom=sqrt((N2*(2*N+3)-sum(ni2.*(2.*ni+3)))/72);
JT=abs(num/denom);
%Compare the JT stats with a standard Normal Distribution
p=1-0.5*erfc(-JT/realsqrt(2)); %p-value
fprintf('JT\t\t\t1-sided p\t\t\tH0\n')
fprintf('%0.4f\t\t %0.4f\t\t\t',JT,p)
if p >= alpha;
   fprintf(' accepted.\n');
else
   fprintf(' rejected.\n');
end    
