clc; 
clf; 
clear all; 
close all; 
n2=input('enter the ending point:n_end='); 
sample_point=0:1:n2; 
amplitude=sample_point; 
subplot(121); 
stem(sample_point,amplitude); 
grid on; 
xlabel('Sample points'); 
ylabel('Amplitude'); 
title('Ramp function'); 
subplot(122); 
plot(sample_point,amplitude); 
grid on; 
xlabel('Sample points'); 
ylabel('Amplitude'); 
title('Ramp function'); 