function LPIIBNRA()
seed=1;
CV=5;
cross_validation(seed,CV);
end

function cross_validation(seed,CV)
    load data.mat;
    rand('seed', seed);
    crossval_idx = crossvalind('Kfold',inter(:),CV);
    final=zeros(990,27);
    for fold=1:5
        fprintf('Round %d validation\n',fold);
        test_idx  = find(crossval_idx==fold);
        matrix = inter;
        matrix(test_idx) = 0;
        rna_sim1=interaction_similarity(matrix,'rna');
        pro_sim1=interaction_similarity(matrix,'pro');
        rna_sim=get_sim(rna_sim1,rna_exp_sim,'rna');
        pro_sim=get_sim(pro_sim1,pro_inter_sim,'pro');
        y_train=recontruction(matrix,rna_sim,pro_sim);
        F = resource_allocate(y_train);
        final(test_idx)=F(test_idx);
    end
    
    [~,~,~,AUC] = perfcurve(inter(:),final(:),1);
    [Recall,Precision,~,AUPR] = perfcurve(inter(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec'); 
    [m,~]=size(Recall);
    F1=[];
    for i=1:m
        F1(i)=(2*Recall(i)*Precision(i))/(Recall(i)+Precision(i));
    end
    [f1,index]=max(F1);
    fprintf('AUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision: %.3f  F1: %.3f \n', AUC, AUPR, Recall(index),Precision(index),f1)
end

function sim=get_sim(sim1,sim2,type)
    if type=='rna'
        sim=(sim1+sim2)/2;
    else
        sim=zeros(27,27);
        for i=1:27
            for j=1:27
                if sim1(i,j)~=0
                    sim(i,j)=(sim1(i,j)+sim2(i,j))/2;
                else
                    sim(i,j)=sim2(i,j);
                end
            end
        end
    end
end

function result = interaction_similarity(inter2,type)
    [rna_num, pro_num] = size(inter2);
    total = sum(sum(inter2));
    if type == 'rna'
        result=zeros(rna_num,rna_num); 
        gama=rna_num/total;
        for i=1:rna_num
            for j=i+1:rna_num
                dis=pdist2(inter2(i,:),inter2(j,:),'squaredeuclidean');
                result(i,j)=exp(-gama*sum(dis));
                result(j,i)=exp(-gama*sum(dis));
            end
        end
    else
        result = zeros(pro_num,pro_num); 
        gama = pro_num/total;
        for i=1:pro_num
            for j=i+1:pro_num
                dis=pdist2(inter2(:,i)',inter2(:,j)','squaredeuclidean');
                result(i,j)=exp(-gama*sum(dis));
                result(j,i)=exp(-gama*sum(dis));
            end
        end
    end
end

function matrix=recontruction(interaction,rna_sim,pro_sim)
    [m,n]=size(interaction);
    result1=zeros(990,27);
    result2=zeros(990,27);
    for i=1:m
        for j=1:m
            result1(i,:)=result1(i,:)+rna_sim(i,j)*interaction(i,:);
        end
        result1(i,:)=result1(i,:)/sum(rna_sim(i,:));
    end
    for i=1:n
        for j=1:n
            result2(:,i)=result2(:,i)+pro_sim(i,j)*interaction(:,i);
        end
        result2(:,i)=result2(:,i)/sum(pro_sim(i,:));
    end
    matrix=(result1+result2)/2;
end

function score_matrix=resource_allocate(interaction_matrix)
    [~,sideffect_num]=size(interaction_matrix);
    resource_allocate_matrix=zeros(sideffect_num,sideffect_num);
    drug_degree=sum(interaction_matrix,2);
    sideffect_degree=sum(interaction_matrix,1);
    for i=1:sideffect_num
        for j=1:sideffect_num
            z=0;
           set1=find(interaction_matrix(:,i)==1);
           set2=find(interaction_matrix(:,j)==1);
           set = intersect(set1,set2);  
           if  ~isempty(set)
              num=size(set,1);
              for p=1:num
                if drug_degree(p)~=0
                  z=z+interaction_matrix(set(p,1),i)*interaction_matrix(set(p,1),j)/drug_degree(set(p,1));
                end
              end
           end

            if sideffect_degree(j)~=0
                resource_allocate_matrix(i,j)=z/sideffect_degree(j);
            end
        end
    end
    score_matrix=((resource_allocate_matrix-0.7*resource_allocate_matrix*resource_allocate_matrix)*interaction_matrix')';
end