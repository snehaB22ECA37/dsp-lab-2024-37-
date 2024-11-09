%{ 
            VERIFICATION OF THE SAMPLING THEOREM WITH A COSINE WAVE UNDER 
                 (I)   UNDER SAMPLING   (Fs < 2*Fm)
                 (II)  AT NYQUIST RATE  (Fs = 2*Fm)
                 (III) AT OVER SAMPLING (Fs >> 2*Fm)
%}

clc;
clf;
close all;

% User Inputs
Fm = input('Enter analog frequency (Hz): ');  
Am = input('Enter amplitude of the cosine wave: ');  
cy = input('Enter number of cycles: '); 

% Parameters
t_final = (1/Fm) * cy;  
sampling_resolution = 0.00002;  

% Time vector for continuous signal
t = 0:sampling_resolution:t_final;  

% Define analog signal
x_t = Am * cos(2*pi*Fm*t); 

% Undersampling condition (fs1 < 2 * Fm)
fs1 = 1.5 * Fm;  
n1 = 0:1/fs1:t_final;  
x_n1 = Am * cos(2*pi*Fm*n1);  

% Nyquist sampling condition (fs2 = 2 * Fm)
fs2 = 2 * Fm;  
n2 = 0:1/fs2:t_final;  
x_n2 = Am * cos(2*pi*Fm*n2); 

% Oversampling condition (fs3 >> 2 * Fm)
fs3 = 20 * Fm;  
n3 = 0:1/fs3:t_final;  
x_n3 = Am * cos(2*pi*Fm*n3);  

% Plotting
figure;

% Plot Input Cosine Wave
subplot(2,2,1);
plot(t, x_t, 'b');
title('Input Cosine Wave');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
ylim([-1.2*Am, 1.2*Am]); % Adjust y-axis limits

% Undersampled Plot
subplot(2,2,2);
plot(t, x_t, 'b');  
hold on;
plot(n1, x_n1, 'r*-');  
title('Undersampled Plot (fs < 2 * fm)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Analog Signal', 'Undersampled Signal', 'FontSize', 6.5); 
grid on;
ylim([-1.2*Am, 1.2*Am]); % Adjust y-axis limits

% Nyquist Plot
subplot(2,2,3);
plot(t, x_t, 'b');
hold on;
plot(n2, x_n2, 'r*-');  
title('Nyquist Plot (fs = 2 * fm)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Analog Signal', 'Nyquist Sampled Signal', 'FontSize', 6.5); 
grid on;
ylim([-1.2*Am, 1.2*Am]); % Adjust y-axis limits 

% Oversampled Plot
subplot(2,2,4);
plot(t, x_t, 'b');  
hold on;
plot(n3, x_n3, 'r*-');  
title('Oversampled Plot (fs >> 2 * fm)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Analog Signal', 'Oversampled Signal', 'FontSize', 6.5); 
grid on;
ylim([-1.2*Am, 1.2*Am]); % Adjust y-axis limits  
