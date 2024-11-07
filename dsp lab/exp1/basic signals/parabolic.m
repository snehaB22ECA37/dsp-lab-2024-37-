clc; 
clf;  
clear all; 
close all; 
n1=input("enter the starting point:n_start="); 
n2=input("enter the ending point:n_end="); 
sample_point=-n1:1:n2; 
amp_1=zeros(1,n1+1); 
for (i=1:n2) 
amp_2(i)=i^2/2; 
end 
Amplitude=[amp_1 amp_2]; 
subplot(121); 
stem(sample_point,Amplitude); 
grid on; 
xlabel('sample points'); 
ylabel('Amplitude'); 
title('discrete'); 
subplot(122); 
plot(sample_point,Amplitude); 
grid on; 
xlabel('sample points'); 
ylabel('Amplitude'); 
title('continuous');