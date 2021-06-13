clear all, close all, clc

%{

HOMEWORK 1

%}

%% Question 1.1    
close all, clear all, clc
format long;

a_1 = single(8888888);
b_1 = single(-8888887);
c_1 = single(0.3333341);

D_1 = a_1 + b_1 + c_1;

E_1 = a_1 + c_1 + b_1;

Error_1 = E_1/D_1;

Error_1_procent = 1 - (E_1/D_1);

Relative_Error_1 = (0.5*10^(-7))/Error_1

Text_1 = ["Question 1"; 'D = ', num2str(D_1); 'E = ', num2str(E_1);'Error = ', num2str(Error_1); 'Error in procent = ', num2str(100*Error_1_procent), '%'];

disp(Text_1)


% Question 1.2

a_2 = 888888888888888;
b_2 = -888888888888887;
c_2 = 0.3333341;

D_2 = a_2 + b_2 + c_2;

E_2 = a_2 + c_2 + b_2; 

Error_2 = D_2/E_2;

Error_2_procent = 1 - (E_2/D_2);

Relative_Error_2 = (0.5*10^(-14))/Error_2   %

Text_2 = ["Question 2"; 'D = ', num2str(D_2); 'E = ', num2str(E_2);'Error = ', num2str(Error_2); 'Error in procent = ', num2str(100*Error_2_procent), '%'];
disp(Text_2)

% Question 1.3

% Machine epsilon single and double precison

m_eps_single = (2^(-23))

m_eps_double = (2^(-52))
% Epsilon1_1=
% 
% Epsil



%% Question 2.1
close all, clear all, clc
format long;

%Forward and central differencing
clear all, close all, clc
format long

x = 0.5;        %in radians
fx = cos(x);    %Function
f_x_prim = -sin(0.5);
f_x_forward = [];
f_x_forward_error = [];
k = 0;
X = [];

% Forward differencing

for A = 0:0.0001:0.1      %stepsize chosen was 0.0001, this effects the precison
    k = k + 1;
    dx = A;      %stepsize
f_x_forward = [f_x_forward, (cos(x+dx) - cos(x)) / (dx)];  

f_x_forward_error = [f_x_forward_error, abs((f_x_forward(k) + sin(x)))];

    X = [X, A];    %used to store the value for dx for plotting
end

% Central differencing
x = 0.5;
t = 0;
X_central = [];
f_x_central_error = [];
f_x_central = [];

for A = 0:0.0001:0.1       %stepsize chosen was 0.0001, this effects the precison
    t = t + 1;
    dx = A;
f_x_central = [f_x_central, (cos(x+dx)-cos(x-dx))/(2*dx)];  

f_x_central_error = [f_x_central_error, abs((f_x_central(t) + sin(x)))];

    X_central = [X_central, A];    %used to store the value for dx for plotting
end

figure('Name', 'Forward and Central error', 'NumberTitle', 'off')
plot(X, f_x_forward_error, X_central,f_x_central_error)
legend('Forward','Central','Location','southeast');
ylabel({'Difference between error and real'});
xlabel({'dx'});
title({'Difference between error and real value for different dx'});




%% 2.2 and 2.3 Trunction error/Actual error
clear all, close all, clc
format long

dx = 0; %stepsize
x = 0.5;        %in radians
fx = cos(x);    %Function
f_x_prim = -sin(0.5);
f_x_forward = [];
Truncation_forward = [];
Truncation_central = [];
k = 0;
X= [];
Total_error_forward = [];
Total_error_central = [];

% Forward and central differencing truncation error

%abs(cos(x)*(1/2)*dx)
%abs(-(1/6)*sin(x)*(dx^2))

for A = logspace(-20, -1, 1000)  %stepsize chosen was 0.0001, this effects the precison
    k = k + 1;
    dx = A;

    Truncation_forward = [Truncation_forward, abs((-cos(x))*(1/2)*dx)];    %Truncation error with taylor serie of the second degree

    Truncation_central = [Truncation_central, abs((-1/6)*sin(x)*(dx^2))];
 
    Total_error_forward = [Total_error_forward, abs(((cos(x+dx) - cos(x)) / (dx)) + sin(x))];
 
    Total_error_central = [Total_error_central,  abs(((cos(x+dx)-cos(x-dx))/(2*dx)) + sin(x))]; 

    X = [X, A];    %used to store the value for dx for plotting
    
end
figure('Name', 'Forward and Central truncation error', 'NumberTitle', 'off');
loglog(X,Truncation_forward,X,Truncation_central)
ylabel({'Truncation Error'});
xlabel({'dx'});
title({'Truncation error depending on step size dx, central and forward difference'});
legend('Forward','Central','Location','southeast');
hold on;

figure ('Name', 'Forward and Central Actual error', 'NumberTitle', 'off');
loglog(X, Total_error_forward, X, Total_error_central);
ylabel({'Actual Error'});
xlabel({'dx'});
title({'Actual error depending on step size dx, central and forward difference'});
legend('Forward','Central','Location','northwest');

