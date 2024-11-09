clc;
clf;
close all;

% Define the input sequences
x = input('Enter the 1st sequence x[n] as a vector (e.g., [1, 2, 3, 4]): ');
h = input('Enter the 2nd sequence h[n] as a vector (e.g., [1, 2, 3, 4]): ');

% Calculate the length of the sequences
N = length(x);

% Zero-padding sequences to length N (if necessary)
x_padded = [x, zeros(1, N - length(x))];
h_padded = [h, zeros(1, N - length(h))];

% Initialize the result
y_circular = zeros(1, N);

% Circular convolution using concentric circles method
for n = 0:N-1
    y_sum = 0;
    for m = 0:N-1
        % Calculate the circular index
        circular_index = mod(n - m, N);
        if circular_index < 0
            circular_index = circular_index + N;
        end
        % Sum the product
        y_sum = y_sum + x_padded(m + 1) * h_padded(circular_index + 1);
    end
    y_circular(n + 1) = y_sum;
end

% Compute circular convolution using MATLAB's built-in function
y_builtin = cconv(x, h, N);

% Display the results
disp('Circular Convolution result using concentric circles method:');
disp(y_circular);
disp('Circular Convolution result using built-in cconv() function:');
disp(y_builtin);

% Plot the sequences and results
figure;

% Plot input sequence x[n]
subplot(3,2,1);
stem(x, 'filled');
title('Input Sequence x[n]');
xlabel('n');
ylabel('x[n]');
grid on;

% Plot input sequence h[n]
subplot(3,2,2);
stem(h, 'filled');
title('Input Sequence h[n]');
xlabel('n');
ylabel('h[n]');
grid on;

% Plot padded x sequence
subplot(3,2,3);
stem(x_padded, 'filled');
title('Padded x[n]');
xlabel('n');
ylabel('x_padded[n]');
grid on;

% Plot padded h sequence
subplot(3,2,4);
stem(h_padded, 'filled');
title('Padded h[n]');
xlabel('n');
ylabel('h_padded[n]');
grid on;

% Plot circular convolution result (concentric circles method)
subplot(3,2,5);
stem(y_circular, 'filled');
title('Circular Convolution y[n] (Concentric Circles)');
xlabel('n');
ylabel('y[n]');
grid on;

% Plot circular convolution result (built-in function)
subplot(3,2,6);
stem(y_builtin, 'filled');
title('Circular Convolution y[n] (Built-in cconv)');
xlabel('n');
ylabel('y[n]');
grid on;
