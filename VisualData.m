clc; %clear the command window
close all; %close all figures
clear %erase all existing variables

%read data from file output.txt
%will make a array full of float values
fileID = fopen('output.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);
fclose(fileID);

%get sixe of file
s = size(data,1);

%extract important constants from file (appended to end)
N = data(s - 3);
T = data(s - 2);
dt = data(s - 1);
dy = data(s);




m = zeros(N,N);
index = 1;

d = cell(N, 1);
for (k = 1:T)
    for (i = 1:N)
        for (j = 1:N)
            m(i,j) = data(index);
            index = index + 1;
        end
    end
    
    d{k} = m;
end


start_phi = d{1,1};

plus_start = 0;
minus_start = 0;

for (i = 1: numel(start_phi))
   element = start_phi(i);
   
   if element < 0
       minus_start = minus_start + 1;
   else 
       plus_start = plus_start + 1;
   end
end


% final_phi = d{T,1};
% 
% plus_end = 0;
% minus_end = 0;
% 
% for (i = 1: numel(final_phi))
%    element = final_phi(i);
%    
%    if element < 0
%        minus_end = minus_end + 1;
%    else 
%        plus_end = plus_end + 1;
%    end
% end
% 
% plus_end - plus_start
% minus_end - minus_start

  for(k = 1:T)
      figure 
      imshow(d{k,1},'InitialMagnification', 200, 'DisplayRange',[-0.005 0.005]);
  end

% for(k = 1:T)
%    figure
%    fname = sprintf('myfile%d.jpg', k);
%    s = imshow(d{k,1},'InitialMagnification', 100, 'DisplayRange',[-0.005 0.005]);
%    t = sprintf('dx=%0.2d, dt=%d, N=%d, T=%d', dx, 5, N, T);
%    title(t)
%    saveas(s,fname);
% 
% end
% 

