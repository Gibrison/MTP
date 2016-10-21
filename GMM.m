function [u,cov_matrix,w] = GMM(X,N,r)
%function for forminng GMM, returns u = final means of clusters,
%cov_matrix = covariance matrix of each cluster in page format
%inputs; X training data, each column is a feature vector
%        N number of cluster required for each GMM / class
%        r the error for binary split in VQ

thresh = 0.5e-010;
ITER = 200;

[m,l] = size(X);
[u,L] = vector_quantizer(X,N,r);

%parameter initialization
for i = 1:N
    ind = find(L==i);
    
    
    
    w(i) = length(ind)/l;        %weight initialization
    for j =1:length(ind)
        matrx{1,i}(:,j) = X(:,ind(j)); %temp variable storing corresponding FV's for a label
    end  
%     cov_matrix(:,:,i) =  vector_cov(matrx{1,i}).*eye(39);
     cov_matrix(:,:,i) =  cov((matrx{1,i})').*eye(m);

    det_cov = det(cov_matrix(:,:,i));
    if det_cov == 0                %singularity check
        str1 =  sprintf('singular matrix for i = %d ',i)
        break;
    end
    %likely hood probabilty calcuation
    P_likelyhood(i,:) = mvnpdf(X',u(:,i)',cov_matrix(:,:,i)); %calculate multivariate distrubution for P
end

% iteration loop
for iter=1:ITER
    
    for j = 1: l  %apriory probability computation
       P_apriori(:,j) = w'.*P_likelyhood(:,j)/sum(w'.*P_likelyhood(:,j));
    end

    w = sum(P_apriori,2)'/l;  %weight update
    temp_u = X*P_apriori';
    temp_var = X.^2*P_apriori';
%     P_square = sum(X.^2);

    for j = 1:N             %mean and covariance update
        u(:,j) = temp_u(:,j)/(w(j)*l);
        cov_matrix(:,:,j) = diag(temp_var(:,j)/(w(j)*l)-u(:,j).^2); 
    end

    P_likelyhood_temp = P_likelyhood;

    for i = 1:N   %likely hood probabilty calcuation update
        det_cov = det(cov_matrix(:,:,i));

        if det_cov == 0                %singularity check
            str2 =  sprintf('singular matrix for i = %d , iter = %d',i,iter)
            break;
        end

        P_likelyhood(i,:) = mvnpdf(X',u(:,i)',cov_matrix(:,:,i)); %calculate multivariate distrubution for P
    end

    if abs(sum(sum(P_likelyhood))-sum(sum(P_likelyhood_temp))) <= thresh
        break
    end
end
str_final =  sprintf('total iterations for convergence = %d',iter)



%%%% my own
% 
% for i=1:l
%     [lol gg(i)]=max(P_apriori(:,i));
% end
% 
% figure;
% for i=1:N
%     
% T=find(gg==i);
%     m(i)=i;
%   D=X(:,T);
% UUU(:,i)=mean(D');
% p=scatter(D(1,:),D(2,:));title('GMM cluster');
% pause(0.5)
% hold on
% %delete(p);
% end
% 
% 
