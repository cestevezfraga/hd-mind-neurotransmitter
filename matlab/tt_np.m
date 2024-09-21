function out=tt_np(data,design)

    % subroutine called by ttest2_np_fdr.m
    % It should be improved by testing for NaN and testing for division by zero at the end
    %
    % Enrico Glerean http://www.glerean.com 
    
    nt=find(design==1);
    as=find(design==2);
    difference=mean(data(:,nt),2) - mean(data(:,as),2);
    s2x=var(data(:,nt),[],2);
    s2y=var(data(:,as),[],2);
    se=sqrt(s2x/length(nt) + s2y/length(as));
    out=difference./se;
