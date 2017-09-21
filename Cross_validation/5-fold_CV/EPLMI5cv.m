function position=EBLMA5cv(x)
%KATZHMDA5ccv(1,1,0.01,2)
%predict disease-related microbe based on KATZHMDA  in the term of 5-fold cross validation
%A: Binary relations between disease and microbe, 1st column:disease, 2nd column:microbe
load('..\..\Data\Expression_profile_data\Known_lncRNA_miRNA_association.mat');
load('..\..\Data\Expression_profile_data\invalid_lnc_expression.mat')
load('..\..\Data\Expression_profile_data\invalid_mi_expression.mat')
load('..\..\Data\Expression_profile_data\invalid_association_expression.mat')
load('..\..\Data\Expression_profile_data\lnc_expression_similarity_matrix.mat')
load('..\..\Data\Expression_profile_data\mi_expression_similarity_matrix.mat')
invalid_lnc_expression=invalid_lnc_expression;
invalid_mi_expression=invalid_mi_expression;
A=Known_lncRNA_miRNA_association;
nl=max(A(:,1));
nm=max(A(:,2));
a=0.5;
b=0.5;
% nl:the number of lncRNAs
% nm:the number of miRNAs
% pp:the number of known lncRNA-miRNA associations
lncrna_similarity=lnc_expression_similarity;
mirna_similarity=mi_expression_similarity;
mean_lnc=mean(lncrna_similarity(:));
mean_mi=mean(mirna_similarity(:));

pp=size(A,1);
original_interaction=zeros(nl,nm);
for i=1:pp
    original_interaction(A(i,1),A(i,2))=1;
end

V = [];
for i=1:size(invalid_association_expression,1)
    V = [V;invalid_association_expression(i,1)];
end
A(V,:)=[]; %remove 2457 associations
pp=size(A,1);
%interaction: adjacency matrix for the lncRNA-miRNA association network
%interaction(i,j)=1 means microbe j is related to disease i
interaction=zeros(nl,nm);
for i=1:pp
    interaction(A(i,1),A(i,2))=1;
end
save interaction interaction;

lncrna_similarity(invalid_lnc_expression,:) = mean_lnc;
lncrna_similarity(:,invalid_lnc_expression) = mean_lnc;
mirna_similarity(invalid_mi_expression,:) = mean_mi;
mirna_similarity(:,invalid_mi_expression) = mean_mi;

lncrna_similarity1=mapminmax(lncrna_similarity,0,1);
mirna_similarity1=mapminmax(mirna_similarity,0,1);

%implement 5-fold cross validation

T=1;
for ccv=1:5
    ccv
    
    load interaction interaction;
    if ccv<5
        AA=A(x((ccv-1)*floor(pp/5)+1:floor(pp/5)*ccv),:);
        % obtain training sample
        for i=1:floor(pp/5)
            interaction(AA(i,1),AA(i,2))=0;
        end
    else
        AA=A(x((ccv-1)*floor(pp/5)+1:pp),:);
        % obtain training sample
        for i=1:pp-floor(pp/5)*4
            interaction(AA(i,1),AA(i,2))=0;
        end
    end
    
    F_1=zeros(nl,nm);
    F_2=zeros(nm,nl);
    SL=lncrna_similarity1*interaction;
    SM=interaction*mirna_similarity1;
    
    % avoid denominator to be 0
    SL(SL==0)=0.00000000001;
    SM(SM==0)=0.00000000001;
    m=nm;
    n=nl;
    parfor i=1:nl
        ff(i,:)=a*sum((repmat(SL(i,:),n,1).*SL)./repmat(sum(SL),n,1),2)'+(1-a)*sum((repmat(SM(i,:),n,1).*SM)./repmat(sum(SM),n,1),2)';
        F_1(i,:)=b*sum((repmat(ff(i,:),m,1).*SL')./repmat(sum(SL,2),1,m)',2)'+(1-b)*sum((repmat(ff(i,:),m,1).*SM')./repmat(sum(SM,2),1,m)',2)';
    end
    
    SL1=SL';
    SM1=SM';
    n1=m;
    m1=n;
    parfor i=1:n1
        ff2(i,:)=a*sum((repmat(SL1(i,:),n1,1).*SL1)./repmat(sum(SL1),n1,1),2)'+(1-a)*sum((repmat(SM1(i,:),n1,1).*SM1)./repmat(sum(SM1),n1,1),2)';
        F_2(i,:)=b*sum((repmat(ff2(i,:),m1,1).*SL1')./repmat(sum(SL1,2),1,m1)',2)'+(1-b)*sum((repmat(ff2(i,:),m1,1).*SM1')./repmat(sum(SM1,2),1,m1)',2)';
    end
    
    F=(F_1+F_2')/2;
    
    num_test=size(AA,1);
    % obtain the score of tested lncRNA-miRNA interaction
    for i=1:num_test
        finalscore(i,1)=F(AA(i,1),AA(i,2));
    end
    
    % make the score of seed lncRNA-miRNA interactions as zero
    for i=1:nl
        for j=1:nm
            if interaction(i,j)==1
                F(i,j)=-99999999999999999999;
            end
        end
    end
    
    for qq=1:num_test
        % obtain the position of tested disease-microbe interaction as variable position(1,ccv),
        ll1=size(find(F>=finalscore(qq)),1);
        ll2=size(find(F>finalscore(qq)),1);
        position(1,T)=ll2+1+(ll1-ll2-1)/2;
        T=T+1;
    end
    
end
save('position.mat','position');
end