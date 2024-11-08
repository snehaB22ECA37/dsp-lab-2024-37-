%{
                    %% EXPERIMENT NUMBER 6
        TITLE : VISUALIZATION OF DFT AND IDFT
%}

clc;
clf;
close all;

N = input('Enter the value for the N point DFT to be found: ');
j = sqrt(-1);

% Input sequence x[n]
x_n = input('Enter the 1st sequence x[n] as a vector (e.g., [1, 2, 3, 4]): ');
L = length(x_n);

% Check if input length is greater than N
if L > N
    disp('Enter a sequence with fewer than N elements.');
    return;
elseif L < N
    % Zero-pad the sequence if it has less than N elements
    x_n = [x_n, zeros(1, N-L)];
end

% Initialize DFT result array
X_k = zeros(1, N);

% Compute the DFT using the formula
for k = 0:N-1
    for n = 0:N-1
        X_k(k+1) = X_k(k+1) + x_n(n+1) * exp((-j*2*pi*k*n)/N);
    end
end

% Verify using MATLAB built-in function
y_builtin = fft(x_n);
disp('The DFT of the sequence x_n using MATLAB built-in fft is:');
disp(y_builtin);

% Display the manually computed DFT sequence
disp('The DFT sequence (manual computation) is:');
disp(X_k);

% Compute and display the magnitude and phase spectra (manual computation)
MagnitudeX_k = abs(X_k);
PhaseX_k = angle(X_k);

disp('The Magnitude sequence (manual computation) is:');
disp(MagnitudeX_k);

disp('The Phase sequence (manual computation) is:');
disp(PhaseX_k);

% Plot the input sequence and DFT results
t = 0:N-1;

figure;
subplot(2,2,1)
stem(t, x_n, 'filled');
title('Input Sequence x[n]');
xlabel('n');
ylabel('x[n]');
%{
subplot(2,2,2)
stem(t, real(X_k), 'filled');
title('Real Part of DFT X[k]');
xlabel('k');
ylabel('Re(X[k])');

subplot(2,2,3)
stem(t, imag(X_k), 'filled');
title('Imaginary Part of DFT X[k]');
xlabel('k');
ylabel('Im(X[k])');
%}

subplot(2,2,2)
plot(t, MagnitudeX_k);
title('Magnitude Spectrum |X[k]|');
xlabel('k');
ylabel('|X[k]|');


subplot(2,2,3);
plot(t, PhaseX_k);
title('Phase Spectrum of X(k)');
xlabel('k'); 
ylabel('Phase(X[k])');

