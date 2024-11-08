%{
                    %% EXPERIMENT NUMBER 7
        TITLE : VERIFICATION OF PROPERTIES OF DFT
%}

%-----------------------------------------------------------------------%

%{
            (1) VERIFICATION OF PERIODICITY PROPERTY  DFT
%}
%{

clc;
clf;
close all;

% Input the sequence from the user
x = input('Enter the sequence as a vector (e.g., [1, 2, 3, 4]):');

% Calculate the length of the input sequence
N = length(x);
% Extend the sequence to two periods
x = [x x];

% Initialize the output array for DFT calculation
y = zeros(1, 2*N);

% Compute the DFT for the two periods
for k = 1 : 2 * N
    for n = 1 : N
        % Calculate the DFT using the formula
        y(k) = y(k) + exp(-1i * 2 * pi * (k-1) * (n-1)/N) * x(n);
    end
end

% Display the DFT results for the 1st period (1 to N)
disp('DFT of input in 1st period (1 to N): ');
disp((y(1 : N))');

% Display the DFT results for the 2nd period (N+1 to 2N)
disp('DFT of input in 2nd period (N+1 to 2N): ');
disp((y(N+1 : 2*N))');

% Verify the periodicity property of the DFT
if abs(y(1 : N) - y(N+1 : 2*N)) < 10^(-10)
    disp('Periodicity Property of DFT is verified');
else
    disp('Periodicity Property of DFT is not verified');
end

%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%{
             (2) VERIFICATION OF LINEARITY PROPERTY OF DFT
%}
%{

clc;
clf;
close all;
disp('---- 2. Linearity Property -----')

% Input sequences and scalar values
x1 = input('Enter the 1st sequence: ');
x2 = input('Enter the 2nd sequence: ');
a = input('Enter 1st scalar value a: ');
b = input('Enter 2nd scalar value b: ');

% DFT function
function X_k = dft(x)
    N = length(x);
    X_k = zeros(1, N); 
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N);
        end
    end
end

% Perform DFT on x1 and x2
y1 = dft(x1);
y2 = dft(x2);

% Check linearity property
Y1 = dft(a * x1 + b * x2); % DFT of a*x1 + b*x2
Y2 = a * y1 + b * y2;      % a*DFT(x1) + b*DFT(x2)

disp('DFT of a*x1 + b*x2');
disp( Y1);

disp('a*DFT(x1) + b*DFT(x2)');
disp(Y2);
% Compare the results
if max(abs(Y1 - Y2)) < 10^(-10)
    disp('Linearity property of DFT is verified');
else
    disp('Linearity property of DFT is not verified');
end
%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%{
             (3) VERIFICATION OF TIME REVERSAL PROPERTY OF DFT
%}
%{

clc;
clf;
close all;
disp('------3. Time Reversal Property----')

% Input sequence
x1 = input('Enter the sequence: ');

% DFT function to compute the Discrete Fourier Transform
function X_k = dft(x)
    N = length(x);             
    X_k = zeros(1, N);         
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N); 
        end
    end
end

% Length of the input sequence
N = length(x1);

% Compute the DFT of the original sequence x1
y1 = dft(x1);

% Generate the time-reversed version of the sequence
n = 1:N-1;
x2(1) = x1(1);                
x2(n+1) = x1(N-(n-1));       

% Compute the DFT of the time-reversed sequence x2
y2 = dft(x2);

% Calculate the expected DFT of the time-reversed sequence using the property
y(1) = y1(1);                 
y(n+1) = y1(N-(n-1));         

% Display the results
disp('DFT of x(n): ');
disp(y1');

disp('DFT of x(N-n): ');
disp(y2');

disp('DFT of x(N-n) by property: ');
disp(y');

% Verify the time reversal property by comparing the two DFT results
if abs(y - y2) < 10^(-10)
    disp('Time reversal property of DFT is verified');
else
    disp('Time reversal property of DFT is not verified');
end
%}


%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%{
             (4) VERIFICATION OF TIME SHIFTING PROPERTY OF DFT
%}
%{

clc;
close all;
disp('------4. Time Shifting Property----')

% Input sequence and shift value
x1 = input('Enter the sequence x1 as a vector (e.g., [1, 2, 3, 4]): ');
m = input('Enter the shift m: ');

% Length of the sequence
N = length(x1);

% User-defined function to compute DFT of a sequence
function X_k = dft(x)
    N = length(x);             
    X_k = zeros(1, N);         
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N); 
        end
    end
end

% Compute DFT of the original sequence x1
y1 = dft(x1);

% Generate the time-shifted version of the sequence using circshift
x2 = circshift(x1', m);  % Shift the sequence by 'm' positions

% Computing the DFT of the time-shifted sequence x2
y2 = dft(x2);

% Calculating the DFT of the time-shifted sequence 
for k = 1:N
    y(k) = y1(k) * exp((-1i * 2 * pi * (k-1) * m) / N);  
end

% Display the results
disp('DFT of x(n): ');
disp(y1');

disp('DFT of x(n-m) by direct method: ');
disp(y2');

disp('DFT of x(n-m) by property: ');
disp(y');

% Verify the time-shifting property by comparing the two DFT results
if max(abs(y - y2)) < 10^(-10)
    disp('Time shifting property of DFT is verified');
else
    disp('Time shifting property of DFT is not verified');
end

%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%{
             (5) VERIFICATION OF FREQUENCY SHIFTING PROPERTY OF DFT
%}
%{

clc;
clf;
close all;
disp('----- 5. Frequency Shifting Property -----');

% Input sequence and frequency shift value
x1 = input('Enter the sequence as a vector (e.g., [1, 2, 3, 4]): ');
m = input('Enter the frequency shift m: ');

% Length of the input sequence
N = length(x1);

% User-defined function to compute the DFT of a sequence
function X_k = dft(x)
    N = length(x);             % Length of the sequence
    X_k = zeros(1, N);         % Initialize DFT result array
    for k = 0:N-1
        for n = 0:N-1
            % Compute the DFT using the standard formula
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N);
        end
    end
end

% Compute DFT of the original sequence x1
y1 = dft(x1);

% Generate the frequency-shifted version of the sequence
for n = 1:N
    % Multiply the sequence by exp(j2πmn/N) to introduce the frequency shift
    x_fshift(n) = exp(1i * 2 * pi * m * (n - 1) / N) * x1(n);
end

% Compute DFT of the frequency-shifted sequence
y2 = dft(x_fshift);

% Compute the DFT of the original sequence and shift it by m in the frequency domain
y = circshift(y1', m);  % Frequency shifting property (circular shift of DFT)

% Display the results
disp('DFT of x(n): ');
disp(y1');

disp('DFT of e^(j2πmn/N) * x(n) by direct method: ');
disp(y2');

disp('DFT of e^(j2πmn/N) * x(n) by property (shifted): ');
disp(y);

% Verify the frequency shifting property by comparing the two DFT results
if max(abs(y' - y2)) < 10^(-10)
    disp('Frequency shifting property of DFT is verified');
else
    disp('Frequency shifting property of DFT is not verified');
end
%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%{
             (6) VERIFICATION OF CIRCULAR CONVOLUTION PROPERTY OF DFT
%}
%{

clc;
clear all;
close all;
disp('------6. Circular Convolution Property-----');
x1 = input('Enter the 1st sequence : ');
x2 = input('Enter the 2nd sequence : ');
N = length(x1);
M = length(x2);

%user defined function to find the dft
function X_k = dft(x)
    N = length(x);
    X_k = zeros(1, N); 
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N);
        end
    end
end

x1 = [x1, zeros(1,M-N)];
x2 = [x2, zeros(1,N-M)];
x = cconv(x1,x2,max(N,M));
y1 = dft(x1);
y2 = dft(x2);
y = y1 .* y2;
Y = dft(x);
disp('DFT of x1 : ');
disp(y1');
disp('DFT of x2 : ');
disp(y2');
disp('DFT of convolution of x1 and x2 : ');
disp(Y');
disp('DFT of convolution of x1 and x2 by property : ');
disp(y');
if abs(y - Y) < 10^(-10)
disp('Circular Convolution property of DFT is verified');
else
disp('Circular Convolution property of DFT is not verified');
end
%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%


%{
            (7) VERIFICATION OF MULTIPLICATION(MODULATION) PROPERTY OF DFT
%}
%{

clc;
clf;
close all;
disp('------ 7. Multiplication (Modulation) Property -----');

% Input the two sequences
x1 = input('Enter the 1st sequence: ');
x2 = input('Enter the 2nd sequence: ');

% Lengths of the sequences
N = length(x1);
M = length(x2);

% Zero-padding the shorter sequence to match the length of the longer sequence
x1 = [x1, zeros(1, M-N)];
x2 = [x2, zeros(1, N-M)];

% Perform pointwise multiplication of the two sequences
x = x1 .* x2;

% User-defined function to compute the DFT of a sequence
function X_k = dft(x)
    N = length(x);               
    X_k = zeros(1, N);           
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N);
        end
    end
end

% Compute the DFTs of the sequences x1 and x2
y1 = dft(x1);
y2 = dft(x2);

% Circular convolution of the DFTs of x1 and x2 divided by N
y = cconv(y1, y2, max(N, M)) / N;

% Compute the DFT of the pointwise multiplied sequence
Y = dft(x);

% Display the results
disp('DFT of x1: ');
disp(y1');

disp('DFT of x2: ');
disp(y2');

disp('DFT of multiplication of x1 and x2 (direct method): ');
disp(Y');

disp('DFT of multiplication of x1 and x2 by property: ');
disp(y');

% Verify the multiplication (modulation) property of DFT by comparing the results
if max(abs(y - Y)) < 10^(-10)
    disp('Multiplication property of DFT is verified');
else
    disp('Multiplication property of DFT is not verified');
end
%}


%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%


%{
            (8) VERIFICATION OF PARSEVALS THEOREM PROPERTY OF DFT
%}
%{

clc;
clf;
close all;
disp('------ 8. Parseval''s Theorem -----');

% Input the sequence
x = input('Enter the sequence: ');

% Length of the sequence
N = length(x);

% User-defined function to compute the DFT of a sequence
function X_k = dft(x)
    N = length(x);               
    X_k = zeros(1, N);           
    for k = 0:N-1
        for n = 0:N-1
            X_k(k+1) = X_k(k+1) + x(n+1) * exp((-1i * 2 * pi * k * n) / N);
        end
    end
end

% Compute the DFT of the sequence
y = dft(x);

% Compute the sum of the squared magnitudes of the sequence in time domain
X = sum(abs(x) .^ 2);

% Compute the sum of the squared magnitudes of the DFT coefficients in frequency domain
Y = sum(abs(y) .^ 2) / N;

% Display results
disp('DFT of x: ');
disp(y');

disp('Sum of |x|^2 (time domain): ');
disp(X);

disp('1/N * Sum of |DFT(x)|^2 (frequency domain): ');
disp(Y);

% Verify Parseval's Theorem by comparing the two sums
if abs(X - Y) < 10^(-10)
    disp('Parseval''s Theorem of DFT is verified');
else
    disp('Parseval''s Theorem of DFT is not verified');
end
%}

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%




