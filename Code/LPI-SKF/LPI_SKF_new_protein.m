function LPI_SKF_new_protein
seed=0;
CV=5;
cross_validation(seed,CV)
end

function final_result=cross_validation(seed,CV)
load data.mat;
lamuda=2^(-3);
interaction_matrix=inter;
[~,col]=size(interaction_matrix);
rand('state',seed);
CV_matrix=ceil(rand(col,1)*CV);
final=zeros(990,27);
for cv=1:CV
    fprintf('round %d validation\n',cv);
    train_interaction_matrix=interaction_matrix;
    test_index=find(CV_matrix==cv);
    train_interaction_matrix(:,test_index)=0;
    
    K2 = [];
    K2(:,:,1)=pro_blast_sim;
    K2(:,:,2)=pro_ctd_lin;
    K_COM2=SKF({K2(:,:,1),K2(:,:,2)},3,5,0.9);
    score_matrix = LapRLS_pro(K_COM2,train_interaction_matrix, lamuda);
    final(:,test_index)=score_matrix(:,test_index);
end
[~,~,~,AUC] = perfcurve(interaction_matrix(:),final(:),1);
[Recall,Precision,~,AUPR] = perfcurve(interaction_matrix(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec'); 
[m,~]=size(Recall);
F1=[];
for i=1:m
    F1(i)=(2*Recall(i)*Precision(i))/(Recall(i)+Precision(i));
end
[f1,index]=max(F1);
fprintf('AUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision: %.3f  F1: %.3f \n', AUC, AUPR, Recall(index),Precision(index),f1)
end


function [LapA] = LapRLS_pro(W2,inter3, lambda)
[~,num_2] = size(inter3);

S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;
d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*inter3';

LapA=A_2';
end

function [W]=SKF(Wall,K,t,ALPHA)
%This program is recoded by reference follow: 
% Wang, B., Mezlini, A. M., Demir, F., Fiume, M., Tu, Z., Brudno, M., et al. (2014). 
% Similarity network fusion for aggregating data types on a genomic scale. Nature Methods, 11, 333�C337
C = length(Wall);
[m,n]=size(Wall{1});

for i = 1 : C
    newW1{i} = Wall{i}./repmat(sum(Wall{i},1),n,1);
end
sumW1 = zeros(m,n);
for i = 1 : C
    sumW1= sumW1 + newW1{i};
end
for i = 1 : C
    Wall{i} = Wall{i}./repmat(sum(Wall{i},1),n,1);
end
for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end
Wsum = zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
for ITER=1:t
    for i = 1 : C           
         Wall{i}=ALPHA*newW{i}*(Wsum - Wall{i})*newW{i}'/(C-1) + (1-ALPHA)*(sumW1 - newW1{i})/(C-1);
    end   
    Wsum = zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
    end     
end
W = Wsum/C;
w = neighborhood_Com(W,K);
W= W.*w;
end

function newW = FindDominateSet(W,K)
[m,n]=size(W);
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
clear IW2;
clear temp;
end

function similarities_N = neighborhood_Com(similar_m,kk)
similarities_N=zeros(size(similar_m));
mm = size(similar_m,1);
for ii=1:mm	
	for jj=ii:mm
		iu = similar_m(ii,:);
		iu_list = sort(iu,'descend');
		iu_nearest_list_end = iu_list(kk);
		ju = similar_m(:,jj);
		ju_list = sort(ju,'descend');
		ju_nearest_list_end = ju_list(kk);
		if similar_m(ii,jj)>=iu_nearest_list_end & similar_m(ii,jj)>=ju_nearest_list_end
			similarities_N(ii,jj) = 1;
			similarities_N(jj,ii) = 1;
		elseif similar_m(ii,jj)<iu_nearest_list_end & similar_m(ii,jj)<ju_nearest_list_end
			similarities_N(ii,jj) = 0;
			similarities_N(jj,ii) = 0;
		else
			similarities_N(ii,jj) = 0.5;
			similarities_N(jj,ii) = 0.5;
        end
    end
end
end