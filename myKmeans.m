function [u,L] = myKmeans(P,k)

[m,n] = size(P);
count = 1;
for i = 1:n
   L(i) = count;
   count = count+1;
   if count == k
       count = 1;
   end
end

for iter = 1:100
Ltemp = L;

u = zeros(m,k);
for i = 1:k
    ind = find(L==i);                                                       %index values corresponnding to each label
    for j =1:length(ind)
        u(:,i) =( u(:,i)*length(ind) + P(:,ind(j)))/length(ind);                 %mean for each label
    end
end

for i =1:n
    for j = 1:k
        euclid(j)= norm(P(:,i)-u(:,j));                                     %new label assignment
    end
    tempind = find(euclid == min(euclid));
    L(i) = tempind(1);
end

if sum(abs(L-Ltemp))<= 1                                                        %termination condition
    break
end







end        
