function [varargout]=perm_test4(action,varargin)

% PERM_TEST permutation test for values of two groups.
%
% USAGE: [p...]=perm_test4('action',arguments);
% 
% There are two actions: 'init' and 'test'.
% 'init': generates permutation vectors. This vector is used for test.
%          perm=perm_test4('init',no_of_left_group,no_of_right_group,...
%                                no_of_permutation);
% 'test': permutation test, You can set a tolerance of the test; 
%         if the p value is smaller than tolernace, H = 1.
%
%         [results]=perm_test4('test',left_group's_data,right_group's_data,...
%                       permutation_vector [,tolerance]);
%
%         Here, left group's data and right_group's data should have
%         the same number of columns, which is to be compared. 
%         So the data should be [subject_no,nodes]. Also you may test
%         other measures at the same time, If you need to test them,
%         Use the 3D array with the form of [subject no,nodes, measures].
%
% Permutation Test, Written by Cheol E Han on Feb 9, 2011, 
%                   Last Modified by Cheol E Han on Feb 12, 2013 
% Related article: Nichols TE, Holmes AP (2001): 
%       Nonparametric permutation tests for functional neuroimaging: 
%                    A primer with examples. Hum Brain Mapp 15:1?25
% 
% See also ttest2, fdr
%



%% Function body
% at least two arguments are needed
if nargin < 2,
    error('Incorrect call to permutation_test: need ACTION');
end

switch lower(action),
    % .........................................................................
    case 'init'
        if (nargin ~= 4)
            error('Incorrect call to permutation test:Init requires the number of left group, the nubmer of right group, and number of permutations');
        end
        no_left   = varargin{1};
        no_right  = varargin{2};
        no_perm   = varargin{3};
        %permutation vector
        all=no_left+no_right;
        perm=zeros(all,no_perm);
        perm(:,1)=1:all;
        for i=2:no_perm
            perm(:,i)=randsample(all,all);
        end
        str.perm=perm;
        str.no_left=no_left;
        str.no_right=no_right;
        str.no_perm=no_perm;
        str.version='perm_test4';
        varargout(1)={str};
        % .........................................................................
    case 'test' % whether
        if (nargin<4 || nargin>5)
            error('Incorrect call to permutation test:Test requires left group, right group and permutation vector (and tolerance)');
        end
        left=varargin{1};   no_left=size(left,1);
        right=varargin{2};  no_right=size(right,1);
        str=varargin{3};    perm=str.perm;              N=str.no_perm;
        if (nargin==5)            tol=varargin{4};
        else                      tol=0.05;         end
        if no_left+no_right~=str.no_left+str.no_right
            error('Incorrect size of permutation vectors')
        end
        if size(left,2)~=size(right,2) | size(left,3)~=size(right,3)
            error('inconsistent number of columns or the measure dimension between two matrix')
        end
        if isfield(str,'version')==0
            error('wrong version: use permutation vector from perm_test4.')
        elseif numel(strfind(str.version,'perm_test4'))==0
            error('wrong version: use permutation vector from perm_test4.')
        end
        perm_temp_left=1:no_left;
        perm_temp_right=(1+no_left):(no_left+no_right);
        % data = [subjects,node measure] or [subjects, edge weights],
        % you can give multiple measure using 3D array.
        all=[left;right];
        no_edges=size(all,2);   no_measures=size(all,3);
        p_two=NaN*ones(no_edges,no_measures,1);
        p_right=NaN*ones(no_edges,no_measures,1);
        z=zeros(N,no_edges,no_measures);
        nnz_edge=zeros(no_measures,1);
        idx=cell(no_measures,1);
        % compute permutation distributino (faster routine when the matrix is sparse)
        for m=1:no_measures;
            % nonzero edges: how many edges were nonzero for each measure.
            idx{m}=find( sum( double(all(:,:,m))~=0 & ~isnan(all(:,:,m)) ,1)>0 );
            nnz_edge(m)=numel(idx{m});
            for k=1:N
                z(k,idx{m},m)=nanmean(all(perm(perm_temp_left,k),idx{m},m))-nanmean(all(perm(perm_temp_right,k),idx{m},m));
            end
        end
        difference=z(1,:,:);
        for m=1:no_measures;
            for k=1:length(idx{m})
                i=idx{m}(k);
                % two-tail
                p_two(i,m)= nnz( abs(z(:,i,m))>=abs(difference(1,i,m)) )  /N;
                % one-tail
                if min(z(:,i,m))==max(z(:,i,m)); p_right(i,m)=0.5;    continue;  end
                if difference(1,i,m)>0
                    p_right(i,m)=nnz(z(:,i,m)>=difference(1,i,m)) /N;
                else
                    p_right(i,m)=1-nnz(z(:,i,m)<=difference(1,i,m)) /N;
                end
            end
        end
        p_left=1-p_right;
        str.p_two=p_two;
        str.p_left=p_left;
        str.p_right=p_right;
        str.z=z;
        str.diff=difference;
        str.no_edges=no_edges;
        str.nnz_edge=nnz_edge;
        varargout(1)={str};
        % .........................................................................
    otherwise
        error('unknown action');
end
