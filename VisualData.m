clc; %clear the command window
close all; %close all figures
clear %erase all existing variables

%read data from file output.txt
%will make a array full of float values
fileID = fopen('output.txt','r');
formatSpec = '%f';
input = fscanf(fileID,formatSpec);
fclose(fileID);

%get sixe of file
s = size(input,1);

%extract important constants from file (appended to end)
N = input(s - 2);
dx = input(s - 1);
dy = input(s);




m = zeros(N,N);

for (i = 1:N )
    
    
    
    
end

