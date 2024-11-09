clc;
clf;
close all;

x = input('Enter the 1st sequence x[n] as a vector (e.g., [1, 2, 3, 4]): ');
h = input('Enter the 2nd sequence h[n] as a vector (e.g., [1, 2, 3, 4]): ');

% Lengths of sequences
N = length(x);
M = length(h);
L = N + M - 1;

% Initialize the output sequence for manual convolution
y_manual = zeros(1, L);

% Perform linear convolution manually
for n = 1:L
    for k = 1:N
        if (n - k + 1 > 0) && (n - k + 1 <= M)
            y_manual(n) = y_manual(n) + x(k) * h(n - k + 1);
        end
    end
end

% Compute linear convolution using MATLAB's built-in function
y_builtin = conv(x, h);

disp(y_manual);
disp(y_builtin);
% Plot the first input sequence
subplot(2, 2, 1);
stem(x, 'filled');
title('Input Sequence x[n]');
xlabel('n');
ylabel('x[n]');
grid on;
ylim([-max(x), max(x) + 1]); % Adjust y-axis limits

% Plot the second input sequence
subplot(2, 2, 2);
stem(h, 'filled');
title('Input Sequence h[n]');
xlabel('n');
ylabel('h[n]');
grid on;
ylim([-max(h), max(h) + 1]); % Adjust y-axis limits

% Plot the output sequence from manual convolution
subplot(2, 2, 3);
stem(y_manual, 'filled');
title('Output Sequence (Manual Convolution)');
xlabel('n');
ylabel('y[n]');
grid on; % Adjust y-axis limits

% Plot the output sequence from built-in function
subplot(2, 2, 4);
stem(y_builtin, 'filled');
title('Output Sequence (Built-in Convolution)');
xlabel('n');
ylabel('y[n]');
grid on;

% Adjust layout
sgtitle('Linear Convolution: Manual vs Built-in');
