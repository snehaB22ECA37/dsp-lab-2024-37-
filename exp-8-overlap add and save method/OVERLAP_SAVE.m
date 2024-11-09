%{
            (1) OVERLAP SAVE METHOD
%}

clc;
clf;
clear;
close all;

% INPUTS
x = input('Enter 1st Sequence X(n)= ');
h = input('Enter 2nd Sequence H(n)= ');
N = input('Enter length of each block L = ');

% PADDING and SIZE ADJUSTMENTS
length_x = length(x);
length_h = length(h);
m = length_h - 1;  % Length for zero-padding

% OVERLAP-SAVE METHOD SETUP
x = [zeros(1,m) x zeros(1,N)];  
h = [h zeros(1,N - length_h)];          
L = N - length_h + 1;          
k = floor((length(x) - m) / L);
disp('number of block of x(n)') 
disp(k);

p = [];

% Perform the block processing (convolution)
for i = 0:k-1
    y = x(i*L+1 : i*L+N);  
    fprintf('Block %d: ', i+1);
    disp(y);               
    q = circonv(y, h);  
    p = [p, q(length_h:end)];
end

% TRIM THE RESULT TO THE LENGTH lx + lh - 1
trimmed_p = p(1:length_x + length_h - 1);

% BUILT-IN CONVOLUTION FOR COMPARISON
conv_builtin = conv(x(m+1 : end-m), h(1 : length_h));  
conv_builtin = conv_builtin(1:length_x + length_h - 1);
fprintf('\nOutput of built-in conv() function:\n');
disp(conv_builtin);

% Print overlap-save result
fprintf('\nUsing Overlap-Save method (trimmed to lx + lh - 1):\n');
disp(trimmed_p);

% PLOT CONVOLVED SIGNAL
stem(trimmed_p, 'filled');
xlabel('n');
ylabel('Amplitude');
title('convolution using overlap save method');

% FUNCTION FOR CIRCULAR CONVOLUTION
function y = circonv(x, h)
    N = length(x);  
    X = fft(x, N);
    H = fft(h, N);
    Y = X .* H;
    y = ifft(Y, N);
end