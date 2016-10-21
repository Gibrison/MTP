clc;
clear all;
close all;

load data2.mat
fott=300.*ones(500,2);
mmm=[C;C+fott];
figure;
scatter(mmm(:,1),mmm(:,2));
C=mmm;





P=C';
k=10;
N=8;
r=0.02;

[u,L]=myKmeans(P,k);

figure;
scatter(u(1,:),u(2,:));title('K means');


[u1,L1] = vector_quantizer(P,N,r)  ;
figure;
[len1 len2]=size(P);
for i=1:N
    l=find(L1==i);
    DD=P(:,l);
    poo=scatter(DD(1,:),DD(2,:));title('vector');
    pause(0.5)
    hold on
end


lll=size(u1)


figure;
scatter(u1(1,:),u1(2,:));title('vector quantization');



 [u2,cov_matrix,w] = GMM(P,N,r);
 
 
figure;
scatter(u2(1,:),u2(2,:));title('GMM Means');

