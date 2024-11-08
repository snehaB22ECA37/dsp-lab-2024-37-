clc;
clf;
close all;

N = input('Enter the value for the N point IDFT to be found: ');
j = sqrt(-1);

% Input sequence X[k]
X_k = input('Enter the sequence X[k] as a vector (e.g., [1, 2, 3, 4]): ');
L = length(X_k);

% Check if input length is greater than N
if L > N
    disp('Enter a sequence with fewer than N elements.');
    return;
elseif L < N
    % Zero-pad the sequence if it has less than N elements
    X_k = [X_k, zeros(1, N - L)];
end

% Initialize IDFT result array
x_n = zeros(1, N);

% Compute the IDFT using the formula
for n = 0 : N-1
    for k = 0 : N-1
        x_n(n+1) = x_n(n+1) + (X_k(k+1) * exp((j*2*pi*k*n)/N));
    end
    % Normalize the result
    x_n(n+1) = x_n(n+1) / N;
end

% Verify using MATLAB built-in function
y_builtin = ifft(X_k);
disp('The IDFT of the sequence X_k using MATLAB built-in ifft is:');
disp(y_builtin);

% Display the manually computed IDFT sequence
disp('The IDFT sequence (manual computation) is:');
disp(x_n);

% Compute and display the magnitude and phase spectra of X[k]
Magnitude_x_n = abs(x_n);
Phase_x_n = angle(x_n);

disp('The Magnitude spectrum of X[k] is:');
disp(Magnitude_x_n);

disp('The Phase spectrum of X[k] is:');
disp(Phase_x_n);

% Plot the input DFT sequence and IDFT results
t = 0 : N-1;

% Plot DFT sequence X(k) 
subplot(2,2,1)
stem(t, X_k, 'filled');
title('Input Sequence X[k]');
xlabel('k');
ylabel('X[k]');
grid on;

% Plot magnitude spectrum of X(k) 
subplot(2,2,2)
plot(t, Magnitude_x_n);
title('Magnitude Spectrum |X[k]|');
xlabel('k');
ylabel('|X[k]|');
grid on;

% Plot phase spectrum of X(k)
subplot(2,2,3);
plot(t, Phase_x_n);
title('Phase Spectrum of X[k]');
xlabel('k'); 
ylabel('Phase(X[k])');
grid on;