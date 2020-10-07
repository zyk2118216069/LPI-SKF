%This program is recoded by reference follow: 
%Zhang, W., Qu, Q., Zhang, Y., Wang, W., 2018b. The linear neighborhood propagation method for predicting long non-coding RNA¨Cprotein interactions. Neurocomputing 273, 526¨C534. https://doi.org/10.1016/j.neucom.2017.07.065
function otherMethod  %%other method
    warning('off');
    seed = 1;
    cross_validation(seed)
end

function cross_validation(seed)
    rand('state',seed);
    load extracted_interaction.txt;
    load extracted_lncRNA_expression.txt
    feature_matrix = extracted_lncRNA_expression';
    interaction_matrix = extracted_interaction;
    final1=zeros(990,27);
    final2=zeros(990,27);
    final3=zeros(990,27);
    final4=zeros(990,27);
    CV=5;
    crossval_idx = crossvalind('Kfold',interaction_matrix(:),CV);

    for k = 1 : CV 
        fprintf('round %d validation\n',k);
        test_idx  = find(crossval_idx==k);
        y_train = interaction_matrix;
        y_train(test_idx) = 0;
        predict_matrix_rwr = RWR_matrix(y_train, 0.1);
        predict_matrix_cf = CF_matrix(y_train);
        predict_matrix_ra = resource_allocate(y_train);
        predict_matrix_hrwr = heterogeneous_matrix(feature_matrix, y_train', 0.9, 0.9, 0.9)';

        final1(test_idx)=predict_matrix_rwr(test_idx);
        final2(test_idx)=predict_matrix_cf(test_idx);
        final3(test_idx)=predict_matrix_ra(test_idx);
        final4(test_idx)=predict_matrix_hrwr(test_idx)';
    end
    get_result('RWR',interaction_matrix,final1);
    get_result('CF',interaction_matrix,final2);
    get_result('LPBNI',interaction_matrix,final3);
    get_result('LPIHN',interaction_matrix,final4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function get_result(model,interaction_matrix,final)
    [~,~,~,AUC] = perfcurve(interaction_matrix(:),final(:),1);
    [Recall,Precision,~,AUPR] = perfcurve(interaction_matrix(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec'); 
    [m,~]=size(Recall);
    F1=[];
    for i=1:m
        F1(i)=(2*Recall(i)*Precision(i))/(Recall(i)+Precision(i));
    end
    [f1,index]=max(F1);
    fprintf('%s£ºAUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision: %.3f  F1: %.3f \n', model,AUC, AUPR, Recall(index),Precision(index),f1)
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
    score_matrix=(resource_allocate_matrix*interaction_matrix')';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score_matrix = RWR_matrix(interaction_matrix, c)
    similarity_matrix = get_similarity_matrix(interaction_matrix);
    %get transformation matrix, i.e. normalize the similarity matrix
    num = size(similarity_matrix, 2);
    for i = 1 : num
        similarity_matrix(i, i) = 0;
    end
    column_sum_matrix = sum(similarity_matrix);
    sum_diagonal_matrix = pinv(diag(column_sum_matrix));
    transformation_matrix =  similarity_matrix * sum_diagonal_matrix;
    %get initial state, same propability for each
    row_sum_interaction_matrix = pinv(sum(interaction_matrix, 2));
    initial_state_matrix = diag(row_sum_interaction_matrix) * interaction_matrix;
    score_matrix = pinv(eye(num) - c * transformation_matrix) * initial_state_matrix;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function score_matrix = CF_matrix(interaction_matrix)
    protein_similarity_matrix = get_similarity_matrix(interaction_matrix);
    %row-mormalize the protein_similarity_matrix
    row_sum_matrix = sum(protein_similarity_matrix, 2);
    sum_diagonal_matrix = pinv(diag(row_sum_matrix));
    row_normalized_protein_similarity_matrix = sum_diagonal_matrix * protein_similarity_matrix;
    score_matrix = row_normalized_protein_similarity_matrix * interaction_matrix;
end

function similarity_matrix = get_similarity_matrix(interaction_matrix)
    %get intersection matrix
    intersection_matrix = interaction_matrix * interaction_matrix';
    %get denominator matrix
    protein_degree_matrix = sum(interaction_matrix, 2);
    denominator_matrix = sqrt(protein_degree_matrix * protein_degree_matrix');
    %calculate similarity_matrix of protein
    similarity_matrix = intersection_matrix ./ denominator_matrix;
    similarity_matrix(isnan(similarity_matrix)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rna_similarity_matrix = get_ran_similarity_matrix(interaction_matrix)
    rna_similarity_matrix = corrcoef(interaction_matrix);
    rna_similarity_matrix(isnan(rna_similarity_matrix)) = 0;
    rna_similarity_matrix = abs(rna_similarity_matrix);
end

function score_matrix = heterogeneous_matrix(feature_matrix, interaction_matrix, gamma, beta, delta)
    protein_similarity_matrix = get_protein_similarity_matrix(interaction_matrix);
    rna_similarity_matrix = get_ran_similarity_matrix(interaction_matrix);
    row_sum_i = sum(interaction_matrix, 2);
    col_sum_i = sum(interaction_matrix);
    %get protein transformation matrix
    protein_num = size(protein_similarity_matrix);
    for i = 1 : protein_num
        protein_similarity_matrix(i, i) = 0;
    end
    row_sum_matrix_p = sum(protein_similarity_matrix, 2);
    sum_diagonal_matrix_p = pinv(diag(row_sum_matrix_p));
    transformation_matrix_p =  sum_diagonal_matrix_p * protein_similarity_matrix .* (1 - gamma);
    transformation_matrix_p(row_sum_i == 0, :) = transformation_matrix_p(row_sum_i == 0, :) ./ (1 - gamma);
    %get lncRNA transformation matrix
    rna_num = size(rna_similarity_matrix);
    for i = 1 : rna_num
        rna_similarity_matrix(i, i) = 0;
    end
    row_sum_matrix_l = sum(rna_similarity_matrix, 2);
    sum_diagonal_matrix_l = pinv(diag(row_sum_matrix_l));
    transformation_matrix_l = sum_diagonal_matrix_l * rna_similarity_matrix .* (1 - gamma);
    transformation_matrix_l(:, col_sum_i == 0) = transformation_matrix_l(:, col_sum_i == 0) ./ (1 - gamma);
    %get protein-lncRNA transformation matrix and lncRNA-protein transformation matrix
    row_sum_matrix_pl = row_sum_i;
    sum_diagonal_matrix_pl = pinv(diag(row_sum_matrix_pl));
    transformation_matrix_pl = sum_diagonal_matrix_pl * interaction_matrix .* gamma;
    %%%
    row_sum_matrix_lp = col_sum_i;
    sum_diagonal_matrix_lp = pinv(diag(row_sum_matrix_lp));
    transformation_matrix_lp = sum_diagonal_matrix_lp * interaction_matrix' .* gamma;
    transformation_matrix_lp(isnan(transformation_matrix_lp)) = 0;
    transformation_matrix_pl(isnan(transformation_matrix_pl)) = 0;
    %get transformation matrix
    transformation_matrix = [transformation_matrix_p transformation_matrix_pl; transformation_matrix_lp transformation_matrix_l];
    %get initial state
    rna_initial_state = eye(rna_num) * (1 - beta);
    protein_initial_state = interaction_matrix *  sum_diagonal_matrix_lp * beta;
    initial_state_matrix = [protein_initial_state; rna_initial_state];
    %get score matrix
    num = size(transformation_matrix, 1);
    score_matrix = (eye(num) / ((eye(num) - (1 - delta) * transformation_matrix'))) * delta * initial_state_matrix;
    score_matrix = score_matrix([1 : protein_num], :);
end

function protein_similarity_matrix = get_protein_similarity_matrix(interaction_matrix)
    %get intersection matrix
    intersection_matrix = interaction_matrix * interaction_matrix';
    %get denominator matrix
    protein_degree_matrix = sum(interaction_matrix, 2);
    denominator_matrix = sqrt(protein_degree_matrix * protein_degree_matrix');
    %calculate similarity_matrix of protein
    protein_similarity_matrix = intersection_matrix ./ denominator_matrix;
    protein_similarity_matrix(isnan(protein_similarity_matrix)) = 0;
end