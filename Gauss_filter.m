%% gaussian filtering

clc;
clear all;
close all;
load nw_data2
M=data_mat(find(data_mat(:,7)==1 & data_mat(:,6)==1),:);
distt=[];
for i=1:76
   t=find(M(:,5)==i);
   if numel(t)>1
   dist=sqrt((M(t(2):t(end),1)-M(t(1):t(end)-1,1)).^2+(M(t(2):t(end),2)-M(t(1):t(end)-1,2)).^2);
   distt=[distt;dist];
   end

end

x(:,1)=1:0.2:10;
sigma = 5;
size = 30;
%x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

yfilt = filter (gaussFilter,5, x);

yfilt2 = conv (x, gaussFilter, 'same');

subplot(2,2,1)
plot(1:46,x)
subplot(2,2,2); plot(1:46,gaussFilter)
subplot(2,2,3); plot(1:46,yfilt)
subplot(2,2,4); plot(1:46,yfilt2)

figure;
subplot(1,2,1);plot(M(:,1),-M(:,2),'.')
% sig=100;
% m1=exp(-M(:,1).^2/(2*sig));
% m1=m1/sum(m1);
% %m1=filter(m1,100,M(:,1));
% 
% m2=exp(-M(:,2).^2/(2*sig));
% m2=m2/sum(m2);
% %m2=filter(m2,100,M(:,2));
% 
% subplot(1,2,2);plot(m1,m2,'.')
% mn1=sum(M(:,1))/2711;
% var1=sum((M(:,1)-sum(M(:,1))/2711).^2)/2711;
% n1=(M(:,1)-mn1)/var1;
% 
% mn2=sum(M(:,2))/2711;
% var2=sum((M(:,2)-sum(M(:,2))/2711).^2)/2711;
% n2=(M(:,2)-mn2)/var2;
% figure; plot(n1,-n2,'.')

% new normalization

sig=1;
win=1;
i1=-sig*win:1:sig*win;
w(:,1)=exp(-i1.^2/(2*sig^2))/sum(exp(-i1.^2/(2*sig^2)));
m1=M;
for i=1+win*sig:2711-win*sig
   m1(i,1)=sum(M(i-win*sig:i+win*sig,1).*w);
   m1(i,2)=sum(M(i-win*sig:i+win*sig,2).*w);
end
 subplot(1,2,2);plot(m1(:,1),-m1(:,2),'.')