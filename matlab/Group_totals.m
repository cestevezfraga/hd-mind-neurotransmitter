
% pmccolgan - Last edit 08.03.14

%Script calculating group totals, averages, threshold averages, permutation
%testing for graph metrics and correlations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=35; %Subject number
loc=75; %threshold taken as 75%

%Patient group totals
name_pat = {XXXXXXXXX','YYYYYYYYY'};
        for s1 = 1 : length(name_pat);
        str = sprintf('Graph_n_metrics_Subj%s.mat',name_pat{s1});
        load(str) % load mat file into workspace
        Tpat_nCC(:,:,s1) = nCC; 
        Tpat_nlambda(:,:,s1) = nlambda;
        Tpat_sigma(:,:,s1) = nCC./nlambda;
        Tpat_T(:,:,s1) = T;
        Tpat_BC(:,:,s1) = BC;
        Tpat_CC(:,:,s1) = CC;
        Tpat_lambda(:,:,s1) = lambda;
        Tpat_deg(:,:,s1) = deg;
        Tpat_stren(:,:,s1) = stren;
        Tpat_Eloc(:,:,s1) = Eloc;
        end

%Patient group averages        
    av_pat_nCC = mean(Tpat_nCC,3);
    av_pat_nlambda = mean(Tpat_nlambda,3);
    av_pat_sigma = av_pat_nCC./av_pat_nlambda;
    av_pat_T = mean(Tpat_T,3);
    av_pat_BC = mean(Tpat_BC,3);
    av_pat_CC = mean(Tpat_CC,3);
    av_pat_lambda = mean(Tpat_lambda,3);
    av_pat_deg = mean(Tpat_deg,3);
    av_pat_stren = mean(Tpat_stren,3);
    av_pat_Eloc = mean(Tpat_Eloc,3);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Control group totals    
name_cont = {'XXXXXXXXX','YYYYYYYYY'};
    for s2 = 1 : length(name_cont);
        str = sprintf('Graph_n_metrics_Subj%s.mat',name_cont{s2});
        load(str) % load mat file into workspace
        Tcont_nCC(:,:,s2) = nCC; 
        Tcont_nlambda(:,:,s2) = nlambda;
        Tcont_sigma(:,:,s2) = nCC./nlambda;
        Tcont_T(:,:,s2) = T;
        Tcont_BC(:,:,s2) = BC;
        Tcont_CC(:,:,s2) = CC;
        Tcont_lambda(:,:,s2) = lambda;
        Tcont_deg(:,:,s2) = deg;
        Tcont_stren(:,:,s2) = stren;
        Tcont_Eloc(:,:,s2) = Eloc;
    end

%Control group averages
    av_cont_nCC = mean(Tcont_nCC,3);
    av_cont_nlambda = mean(Tcont_nlambda,3);
    av_cont_sigma = av_cont_nCC./av_cont_nlambda;
    av_cont_T = mean(Tcont_T,3);
    av_cont_BC = mean(Tcont_BC,3);
    av_cont_CC = mean(Tcont_CC,3);
    av_cont_lambda = mean(Tcont_lambda,3);
    av_cont_deg = mean(Tcont_deg,3);
    av_cont_stren = mean(Tcont_stren,3);
    av_cont_Eloc = mean(Tcont_Eloc,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%Threshold averages transposed for ease of results reporting
    cont_deg = av_cont_deg;
    cont_deg = transpose(cont_deg);
    cont_stren = av_cont_stren;
    cont_stren = transpose(cont_stren);
    cont_BC = av_cont_BC;
   
    cont_Eloc = av_cont_Eloc;
   
    cont_CC = av_cont_CC;
    
    
    pat_deg = av_pat_deg;
    pat_deg = transpose(pat_deg);
    pat_stren = av_pat_stren;
    pat_stren = transpose(pat_stren);
    pat_BC = av_pat_BC;
    
    pat_Eloc = av_pat_Eloc;
  
    pat_CC = av_pat_CC;
   
results=cat(2,cont_deg,pat_deg,cont_stren,pat_stren,cont_BC,pat_BC,cont_Eloc,pat_Eloc,cont_CC,pat_CC);

save 'group_results'

