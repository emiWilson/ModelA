clc; %clear the command window
close all; %close all figures
clear %erase all existing variables

fileID = fopen('output.txt','r');

formatSpec = '%f';

A = fscanf(fileID,formatSpec);