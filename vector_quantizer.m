function [u,L] = vector_quantizer(P,N,r)  

%P should be matrix of COLUMN feature vectors 
%N is number of means/classes,done using binary split, N  should be power of 2
% r is the constant used for binary splitting of means 

%[u,L] output correspond to u(Nx1) mean vector and L which is corresponding
%label for each poit

m = size(P,1);  %dimension of feature vectors
l = size(P,2);  % number of feature vectors
u(:,1) = sum(P,2)/l; %initial mean
e = r*ones(m,1);

ITER = 200;     %not passed as a parameter, change withing function, as 200 is a high enough iteration value 

for t =1:log2(N)
    
    k = 2^t;
    
    for i =1:k/2
        utemp(:,i) = u(:,i)+e;
        utemp(:,i+k/2) = u(:,i)-e;
    end
    
    u = utemp; 
    
    for i =1:l
        for j = 1:k
            euclid(j)= norm(P(:,i)-u(:,j));                                     %new label assignment
        end
        tempind = find(euclid == min(euclid));
        L(i) = min(tempind);
    end

    %mean optimization start here
    %*************************************************************************************************************%
%     u = zeros(m,k);
    for iter = 1:ITER
        
 
        
        Ltemp = L;
        
        
        
        for i = 1:k
            ind = find(L==i); 
           %index values corresponnding to each label
           TF=isempty(ind);
           if TF==0
               
       
           
           sum_u = 0;
            for j =1:length(ind)
                sum_u =( sum_u*length(ind) + P(:,ind(j)))/length(ind);                 %mean for each label
            end
            u(:,i) = sum_u;
           end
        end
        
        for i =1:l
            for j = 1:k
                euclid(j)= norm(P(:,i)-u(:,j));                                     %new label assignment
            end
            tempind = find(euclid == min(euclid));
            L(i) =min(tempind);
        end

        if sum(abs(L-Ltemp))== 0                                                        %termination condition
            break;
        end


    end
    

    
    
  
    
end


