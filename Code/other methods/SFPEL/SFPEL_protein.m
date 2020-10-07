%This program is recoded by reference follow: 
%Zhang, W., Yue, X., Tang, G., Wu, W., Huang, F., Zhang, X., 2018c. SFPEL-LPI: Sequence-based feature projection ensemble learning for predicting LncRNA-protein interactions. PLoS Comput Biol 14. https://doi.org/10.1371/journal.pcbi.1006616
function SFPEL_protein
seed=0;
CV=5;
cross_validation(seed,CV)
end

function cross_validation(seed,CV)
    load('NP2.0.mat')
    interaction_matrix=InteractionMatrix;
    [~,col]=size(interaction_matrix);
    rand('state',seed);
    CV_matrix=ceil(rand(col,1)*CV);
    final=zeros(990,27);
    for cv=1:CV
        fprintf('round %d validation\n',cv);
        train_interaction_matrix=interaction_matrix;
        train_index=find(CV_matrix~=cv);
        test_index=find(CV_matrix==cv);
        train_interaction_matrix(:,test_index)=0;

        X_feature_num=2;
        train_X_cell{1,1}=PCPseAACFeature_Protein;
        test_X_cell{1,1}=PCPseAACFeature_Protein;
        train_X_cell{1,1}(test_index,:)=0;
        test_X_cell{1,1}(train_index,:)=0;

        train_X_cell{2,1}=SCPseAACFeature_Protein;
        test_X_cell{2,1}=SCPseAACFeature_Protein;
        train_X_cell{2,1}(test_index,:)=0;
        test_X_cell{2,1}(train_index,:)=0;

        for j=1:X_feature_num
            train_protein_simi_cell{j,1}=fast_LNS_caculate(train_X_cell{j,1},size(train_X_cell{j,1},1));
        end
        score_matrix=LPI_pred(test_X_cell,train_interaction_matrix',train_protein_simi_cell,test_X_cell); %prediction results for new lncRNAs
        score_matrix=score_matrix';
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
    fprintf('AUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision:: %.3f  F1: %.3f \n', AUC, AUPR, Recall(index),Precision(index),f1)
end

function y_pred=LPI_pred(Xs,train_interaction_matrix,train_lncRNA_simi_cell,test_X_cell)
% initialize the parameters in the objective fuction
mu=0.001;
lam=0.0001;
gam=4;
max_iter=10;
eps=1e-6;
Y=train_interaction_matrix;
simi_num=size(train_lncRNA_simi_cell,1);
for i=1:simi_num
    Ls{i,1}=eye(size(train_lncRNA_simi_cell{i,1},1))-train_lncRNA_simi_cell{i,1};
end
mx=size(Xs,1);
ml=size(Ls,1);
for i=1:mx
    Xs{i,1}=Xs{i,1}';
end
% predict interactions for the lncRNAs in the testing dataset (y_pred)
train_res=lpi(Xs, Ls, Y, mx, ml, mu, lam, gam, max_iter, eps);
Gs=train_res.Gs;
alpha=train_res.alpha;
y_pred=zeros(size(test_X_cell{1,1},1),size(train_interaction_matrix,2));
for i=1:mx
    y_pred=y_pred+alpha(i,1)*test_X_cell{i,1}*Gs{i,1};
end
end

function train_res=lpi(Xs, Ls, Y, mx, ml, mu, lam, gam, max_iter, eps)
% optimize the objective fuction
[n,c]=size(Y);
Gs=cell(mx,1);
e1ds=cell(mx,1);
ds=zeros(mx,1);
for i=1:mx
    ds(i,1)=size(Xs{i,1},1);
    Gs{i,1}=rand(ds(i,1),c);
    e1ds{i,1}=ones(ds(i,1),1);
end
alpha=ones(ml,1)/ml;
Gs_old=cell(mx,1);
As=cell(mx,1);
As_pos =cell(mx,1);
As_neg =cell(mx,1);
Bs =cell(mx,1);
Bs_pos =cell(mx,1);
Bs_neg =cell(mx,1);
L=zeros(n,n);
for i=1:ml
    L=L + alpha(i,1).^gam*Ls{i,1};
end
t=0;
while(t < max_iter)
    t=t + 1;
    Q =Y;
    for i=1:mx
        Gs_old(i,1)= Gs(i,1);
        Q=Q + mu*Xs{i,1}'*Gs{i,1};
    end
    P=inv(L + (1 + mx*mu)*diag(ones(n,1)));
    F_mat=P*Q;
    for i=1:mx
        As{i,1}=Xs{i,1}*(mu*diag(ones(n,1))-mu^2*P')*Xs{i,1}' + lam*(e1ds{i,1}*e1ds{i,1}');
        As_pos{i,1}=(As{i,1} + abs(As{i,1}))/2;
        As_neg{i,1}=(abs(As{i,1}) - As{i,1})/2;
        Bs{i,1}=mu*Xs{i,1}*P*Y;
        for j=1:mx
            if(i == j)
                continue;
            else
                Bs{i,1}=Bs{i,1} + mu^2*Xs{i,1}*P'*Xs{j,1}'*Gs{j,1};
            end
        end
        Bs_pos{i,1}=(Bs{i,1} + abs(Bs{i,1}))/2;
        Bs_neg{i,1}=(abs(Bs{i,1}) - Bs{i,1})/2;
    end
    for i=1:mx
        Gs{i,1}= Gs{i,1}.*sqrt((Bs_pos{i,1} + As_neg{i,1}*Gs{i,1})./(Bs_neg{i,1} + As_pos{i,1}*Gs{i,1}));
    end
    for i= 1:ml
        alpha(i,1)=(1/sum(diag((F_mat)'*Ls{i,1}*F_mat))).^(1/(gam - 1));
    end
    alpha=alpha/sum(alpha);
    L=zeros(n,n);
    for i=1:ml
        L=L + alpha(i,1).^gam*Ls{i,1};
    end
    diff_G=zeros(mx,1);
    for i=1:mx
        diff_G(i,1)=norm(Gs{i,1} - Gs_old{i,1}, 'fro')/norm(Gs_old{i,1},  'fro');
    end
    if(mean(diff_G) < eps)
        break
    end
end

for i=1:mx
    train_res.Gs{i,1}=Gs{i,1};
end
train_res.alpha=alpha;
train_res.predict=F_mat;
end



function W=fast_LNS_caculate(feature_matrix,neighbor_num)
%get the linear neighborhood similarities
iteration_max=50;
mu=0;
X=feature_matrix;

row_num=size(X,1);
distance_matrix=pdist2(X,X,'euclidean');
e=ones(row_num,1);
distance_matrix=distance_matrix+diag(e*inf);
[~, si]=sort(distance_matrix,2,'ascend');
nearst_neighbor_matrix=zeros(row_num,row_num);
index=si(:,1:neighbor_num);
for i=1:row_num
    nearst_neighbor_matrix(i,index(i,:))=1;
end

C=nearst_neighbor_matrix;
rand('state',2);
W=rand(row_num,row_num);
W=(C.*W);
lamda=8*e;
P=X*X'+lamda*e';
for i=1:iteration_max
    Q=(C.*W)*P+mu*(C.*W);
    W=(C.*W).*P./Q;
    W(isnan(W))=0;
end
end