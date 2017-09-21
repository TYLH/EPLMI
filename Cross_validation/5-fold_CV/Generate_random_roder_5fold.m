load('..\..\Data\Expression_profile_data\Known_lncRNA_miRNA_association.mat')
load('..\..\Data\Expression_profile_data\invalid_association_expression.mat')

V = [];
for i=1:size(invalid_association_expression,1)
    V = [V;invalid_association_expression(i,1)];
end
Known_lncRNA_miRNA_association(V,:)=[]; %remove 2457 associations

num=size(Known_lncRNA_miRNA_association,1);

for cv=1:20
    x=randperm(num)';
    Random_order(cv,:)=x;
end

save Random_order.mat Random_order
