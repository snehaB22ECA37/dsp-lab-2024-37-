%{
            (1) OVERLAP ADD METHOD
%}

clc;  
clf;  
clear;  
close all;  

% USER INPUT
x = input('Enter the First Sequence, x: = ');  
h = input('Enter the Second Sequence, h: = ');  

% LENGTHS OF INPUT SEQUENCES
h1 = length(x);  
h2 = length(h);  

% ZERO-PADDING THE FIRST SEQUENCE
x = [x, zeros(1, h2)];  

% CHECK IF THE FIRST SEQUENCE IS LONGER THAN THE SECOND
if (h1 > h2)
    s = [zeros(1, h1 + (2 * h2) - 1)];
    
    % PERFORM OVERLAP-ADD METHOD
    for (i = 1:h2:h1)
        a = x(i:i + h2 - 1);  
        s1 = convdft(a, h);  
        s(i:i + (2 * h2) - 2) = s(i:i + (2 * h2) - 2) + s1;  
    end
end

% overlap add method output
y = s(1:h1 + h2 - 1);  
disp('LC using overlap add method');
disp(y);  

% builtin function output
conv_builtin = conv(x, h);  
conv_builtin = conv_builtin(1:h1 + h2 - 1);
disp('LC using builtin function');
disp(conv_builtin); 

%Plotting the output sequences
subplot(2,1,1)
stem(y,'filled')
xlabel('n--->');
ylabel('Amplitude --->');
title('Convolution using overlap add method')

subplot(2,1,2)
stem(conv_builtin,'filled')
xlabel('n--->');
ylabel('Amplitude --->');
title('Convolution using builtin function')

% FUNCTION FOR DFT-BASED CONVOLUTION
function h = convdft(x, y)
    n1 = length(x);  
    n2 = length(y);    
    x = [x, zeros(1, n2 - 1)];  
    y = [y, zeros(1, n1 - 1)];  
    dx = fft(x);  
    dy = fft(y);  
    mul = dx .* dy;   
    h = ifft(mul);  
end