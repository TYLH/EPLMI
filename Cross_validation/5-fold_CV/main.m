Generate_random_roder_5fold;
for cv=1:20
cv
load('Random_order.mat')
p=EPLMI5cv(Random_order(cv,:));
POSITION(cv,:)=p;
Overallauc(cv)=Position2AUC(p);
end
save Overallauc Overallauc
save POSITION POSITION
a=mean(Overallauc);
b=std(Overallauc);