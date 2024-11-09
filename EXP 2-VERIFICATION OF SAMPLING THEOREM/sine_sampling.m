%{ 
            VERIFICATION OF THE SAMPLING THEOREM WITH A SINE WAVE UNDER 
                 (I)   UNDER SAMPLING   (Fs < 2*Fm)
                 (II)  AT NYQUIST RATE  (Fs = 2*Fm)
                 (III) AT OVER SAMPLING (Fs >> 2*Fm)
%}

clc;
clf;
close all;

% User Inputs
Fm = input('Enter the frequency of the sine wave (Hz): ');  
Am = input('Enter the amplitude of the sine wave: ');       
Cy = input('Enter the number of cycles: ');

% Parameters
t_final = (1/Fm) * Cy;  
sampling_resolution = 0.0002;  

% Time vector for continuous signal
t = 0:sampling_resolution:t_final;  

% Define sine wave signal
x_t = Am * sin(2*pi*Fm*t); 

% Undersampling condition (fs1 < 2*fm)
fs1 = 1.5 * Fm;  
n1 = 0:(1/(fs1)):t_final;  
x_n1 = Am * sin(2*pi*Fm*n1);  

% Nyquist sampling condition (fs1 = 2*fm)
fs2 = 3 * Fm;  
n2 = 0:(1/(fs2)):t_final;  
x_n2 = Am * sin(2*pi*Fm*n2);  

% Oversampling condition (fs1 >> 2*fm)
fs3 = 2000 * Fm;  
n3 = 0:(1/(fs3)):t_final;  
x_n3 = Am * sin(2*pi*Fm*n3);  

% Plotting
figure;

% Plot Input Sine Wave
subplot(2,2,1);
plot(t, x_t, 'b');
title('Input Sine Wave');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
ylim([-1.2*Am, 1.2*Am]);

% Undersampled Plot
subplot(2,2,2);
plot(t, x_t, 'b');
hold on;
plot(n1, x_n1, 'r*-');
title('Undersampled Plot (Fs < 2*Fm)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
ylim([-1.2*Am, 1.2*Am]);

% Nyquist Plot
subplot(2,2,3);
plot(t, x_t, 'b');
hold on;
plot(n2, x_n2, 'r*-');
title('Nyquist Plot (Fs = 3*Fm)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Continuous Signal', 'Nyquist Sampled Signal', 'FontSize', 6.5);
grid on;
ylim([-1.2*Am, 1.2*Am]);

% Oversampled Plot
subplot(2,2,4);
plot(t, x_t, 'b');
hold on;
plot(n3, x_n3, 'r*-');
title('Oversampled Plot (Fs >> 2*Fm)' );
xlabel('Time (s)');
ylabel('Amplitude');
legend('Continuous Signal', 'Oversampled Signal', 'FontSize', 6.5);
grid on;
ylim([-1.2*Am, 1.2*Am]);
