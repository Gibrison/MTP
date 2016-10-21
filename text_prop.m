%% text properties in matlab

clc;
clear all;
close all;

t = text(0.5,0.5,'hi,\n  ');
s = t.FontSize;
t.FontSize = 12;
'my label';

{'first line','second line'};
disp('K_{n}');

sprintf('\color{magenta} text');