clc; 
clf; 
close all; 
clear all; 
n1=input('Enter the starting point of signal: n_start'); 
n2=input('Enter the end point of signal: n_end'); 
sample_points = -n1:1:n2; 
amplitude = [ zeros(1,n1) 1 ones(1,n2)]; 
stem(sample_points, amplitude) 
grid on; 
xlabel('samples'); 
ylabel('amplitude');