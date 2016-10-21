
%% step: feature matrix generation 
%=====================================================================
clc;
close all;
clear all;
load('denoised_M.mat')
%laod data_mat
%laod A
load A_new
%load M1
A( ~any(A,2), : ) = [] ;%rows
D_mat=M;
D_mat(:,4)=A;  % raw data matrix formed

cut=find(D_mat(:,7)==20);
% for 20 writers
Mat_20=D_mat(1:cut(end),:);
% %normalization
[ro,co]=size(Mat_20);
mx=sum(Mat_20(:,1))/ro;
sig_x=sqrt(sum((Mat_20(:,1)-mx).*(Mat_20(:,1))/ro));
Mat_20(:,1)=(Mat_20(:,1)-mx)/sig_x;

my=sum(Mat_20(:,2))/ro;
sig_y=sqrt(sum((Mat_20(:,2)-my).*(Mat_20(:,2))/ro));
Mat_20(:,2)=(Mat_20(:,2)-my)/sig_y;
label(1:680858,1)=Mat_20(1:ro-1,7);


%%==========================================
%% point based features
% velocity , writing direction ,curvature.

for i=1:20
    len=find(Mat_20(:,7)==i);
    %sz=size(len);
    label(len(end))=0;
    time(len(1):len(end)-1,1)=Mat_20(1+len(1):len(end),3)-Mat_20(len(1):len(end)-1,3);
   li(len(1):len(end)-1,1)=sqrt((Mat_20(1+len(1):len(end),1)-Mat_20(len(1):len(end)-1,1)).^2+(Mat_20(1+len(1):len(end),2)-Mat_20(len(1):len(end)-1,2)).^2);
   cos_t(len(1):len(end)-1,1)=Mat_20(1+len(1):len(end),1)-Mat_20(len(1):len(end)-1,1);
   sin_t(len(1):len(end)-1,1)=Mat_20(1+len(1):len(end),2)-Mat_20(len(1):len(end)-1,2);
   cos_phi(len(1):len(end)-2,1)=cos_t(len(1):len(end)-2,1).*cos_t(len(2):len(end)-1,1)+sin_t(len(1):len(end)-2,1).*sin_t(len(2):len(end)-1,1);
   sin_phi(len(1):len(end)-2,1)=cos_t(len(1):len(end)-2,1).*sin_t(len(2):len(end)-1,1)-sin_t(len(1):len(end)-2,1).*cos_t(len(2):len(end)-1,1);
   
end
label=label(1:ro-1);
cos_phi(ro-1,1)=0;
sin_phi(ro-1,1)=0;
%velocity=li./time;
FM=zeros(ro-1,6);
FM(:,1)=li./time;
%zeros=find(li==0);
li(find(li==0))=1;
%cos_teta=cos_t./li;
%sin_teta=sin_t./li;

FM(:,2)=cos_t./li;
FM(:,3)=sin_t./li;
FM(:,4)=cos_phi./(li.*li);
FM(:,5)=sin_phi./(li.*li);
FM(:,6)=label;
FM( ~any(FM,2), : ) = [] ;%rows
