%{

Homework 2

%}

%% Question 1
clear all, close all, clc
format long;

%Given data:
L = 0.4;     %[m] x-direction             
H = 0.2;     %[m] y-direction
k = 54.3;    %[W/m*K]  thermal conductivity
Cp = 460;    %[J/kg*K] specific heat capacity
rho = 7800;  %[kg/m^3] density

%From Homework 1 3 point central dx = dy = 0.0125 m [Direct method]

T_point_1 = 19.8569;
T_point_2 = 20.525;
T_point_3 = 37.3569;

alpha = k/(rho*Cp); %[m^2/s] thermal diffusitivity

dx = 0.0125;
dy = dx;
x1 = (0.1/dx)+1;
x2 = (0.2/dx)+1;
x3 = round((0.3/dx)+1);
y_point = (0.1/dy) +1;

    %[vertical(y), horizontal(x)]
            
%Amounts of points in the Rows/collums
M = (H/dy)+1;       % y
N = (L/dx)+1;       % x 
            
%   A)
%Makes the matrices for both step sizes
T = zeros(M, N);  

%Inputs the boundary condition (Corner pieces not correct = 0)
Top_BC = 20;
Bottom_BC = 10;
Left_BC = 30;
Right_BC = 100;
%BOUNDARY CONDITION APPLIED TO THE MATRIX
T(M,2:N - 1) = Bottom_BC;
T(1,2:N - 1) = Top_BC;
T(2:M -1,1) = Left_BC;
T(2:M -1,N) = Right_BC;

%Gives average values for corner pieces, just for the 2d contur plot
T(1,1) = (Left_BC + Top_BC)/2;
T(1,N) = (Right_BC + Top_BC)/2;
T(M,1) = (Left_BC + Bottom_BC)/2;
T(M,N) = (Right_BC + Bottom_BC)/2;

% EULER METHOD:

tol = 1e-5;     % Given in the question description
error = 1;      % Equals to 1 so it will enter the while loop
iterations = 0;
dt = 4.47;      %Manually tested for lowest iterations required
dt_error = 1;
iterations_new = 10000;

% 
% while dt_error > 0.5
% dt = dt + 3.47 ;            %Mannualy found value for low iterations 4.47
% iterations_before = iterations_new;
    while error > tol
    
    iterations = iterations + 1;
    T_old = T;                 
    Fo_x = (alpha*dt)/(dx^2);
    Fo_y = (alpha*dt)/(dy^2);

      for y = 2:M-1         % In y direction
        for x = 2:N-1     % In x direction
            
            T(y,x) = T(y,x) + Fo_x*(T(y,x-1)+T(y,x+1)-2*T(y,x)) + Fo_y*(T(y-1,x)+T(y+1,x)-2*T(y,x));
            
        end
        
     end

    error = max(max(abs(T-T_old)));
    end
%     iterations_new = iterations;
% dt_error = abs(iterations_new - iterations_before);
% end

X1 = ['At first point ', num2str(T(y_point,x1)), ' Celcius']
X2 = ['At second point ', num2str(T(y_point,x2)), ' Celcius']
X3 = ['At third point ', num2str(T(y_point,x3)), ' Celcius']

%Difference between direct method and time stepping

Diff_X1_Euler = T_point_1 - T(y_point,x1)
Diff_X2_Euler = T_point_2 - T(y_point,x2)
Diff_X3_Euler = T_point_3 - T(y_point,x3)

Iterations_required = iterations

dx2 = dx^2;
dy2 = dy^2;


% max_dt = (1/alpha)*(1/2)*(dx2*dy2)*(1/(dx2+dy2))

max_dt_2 = dx2/(2*alpha)


%% Question 1 RUNGE-KUTTA method  For dx = dy
clear all, close all, clc
format long;

%Given data:
L = 0.4;     %[m] x-direction             
H = 0.2;     %[m] y-direction
k = 54.3;    %[W/m*K]  thermal conductivity
Cp = 460;    %[J/kg*K] specific heat capacity
rho = 7800;  %[kg/m^3] density

%From Homework 1 3 point central dx = dy = 0.0125 m [Direct method]

T_point_1 = 19.8569;
T_point_2 = 20.525;
T_point_3 = 37.3569;

alpha = k/(rho*Cp); %[m^2/s] thermal diffusitivity

dx = 0.0125;
dy = dx;
x1 = (0.1/dx)+1;
x2 = (0.2/dx)+1;
x3 = round((0.3/dx)+1);
y_point = (0.1/dy) +1;

    %[vertical(y), horizontal(x)]
            
%Amounts of points in the Rows/collums
M = (H/dy)+1;       % y
N = (L/dx)+1;       % x 
            
%   A)
%Makes the matrices for both step sizes
T = zeros(M, N);  

%Inputs the boundary condition (Corner pieces not correct = 0)
Top_BC = 20;
Bottom_BC = 10;
Left_BC = 30;
Right_BC = 100;
%BOUNDARY CONDITION APPLIED TO THE MATRIX
T(M,2:N - 1) = Bottom_BC;
T(1,2:N - 1) = Top_BC;
T(2:M -1,1) = Left_BC;
T(2:M -1,N) = Right_BC;

%Gives average values for corner pieces, just for the 2d contur plot
T(1,1) = (Left_BC + Top_BC)/2;
T(1,N) = (Right_BC + Top_BC)/2;
T(M,1) = (Left_BC + Bottom_BC)/2;
T(M,N) = (Right_BC + Bottom_BC)/2;


tol = 1e-5;     % Given in the question description
error = 1;      % Equals to 1 so it will enter the while loop
iterations = 0;
dt = 2.753;             %Manually tested for lowest iterations 2.753
k1 = zeros(M,N);
k2 = zeros(M,N);
k3 = zeros(M,N);
k4 = zeros(M,N);
%Runge kutta 4th order to solve unknown temps in Matrix

while error > tol
    iterations = iterations + 1;
    T_old = T;
    Fo = (alpha*dt)/(dx^2);
    
   
     for y = 2:M-1         % In y direction
        for x = 2:N-1     % In x direction
        
        k1(y,x) = Fo*(T(y-1,x) + T(y+1,x) + T(y,x-1) + T(y,x+1) -4*T(y,x));
        k2(y,x) = Fo*((T(y-1,x) + 0.5*k1(y-1,x)) + (T(y+1,x) + 0.5*k1(y+1,x)) + (T(y,x-1) + 0.5*k1(y,x-1)) + (T(y,x+1) + 0.5*k1(y,x+1)) -(4*T(y,x) + 0.5*k1(y,x)));
        k3(y,x) = Fo*((T(y-1,x) + 0.5*k2(y-1,x)) + (T(y+1,x) + 0.5*k2(y+1,x)) + (T(y,x-1) + 0.5*k2(y,x-1)) + (T(y,x+1) + 0.5*k2(y,x+1)) -(4*T(y,x) + 0.5*k2(y,x)));
        k4(y,x) = Fo*((T(y-1,x) + k3(y-1,x)) + (T(y+1,x) + k3(y+1,x)) + (T(y,x-1) + k3(y,x-1)) + (T(y,x+1) + k3(y,x+1)) -(4*T(y,x) + k3(y,x)));

        T(y,x) = T(y,x) + (1/6)*(k1(y,x) + 2*k2(y,x) + 2*k3(y,x) + k4(y,x));
           
        end
        
     end

    error = max(max(abs(T-T_old)));
end


X1 = ['At first point ', num2str(T(y_point,x1)), ' Celcius']
X2 = ['At second point ', num2str(T(y_point,x2)), ' Celcius']
X3 = ['At third point ', num2str(T(y_point,x3)), ' Celcius']

%Difference between direct method and time stepping

Diff_X1_Euler = T_point_1 - T(y_point,x1)
Diff_X2_Euler = T_point_2 - T(y_point,x2)
Diff_X3_Euler = T_point_3 - T(y_point,x3)

Iterations_required = iterations

contourf(T,20)

%% Question 2 LAX-FREDRICHS METHOD ADVECTION EQUATION
clear all, close all, clc
Points = 100; %in the rod
Length = 1; %m
U = 1;  %ms/s
t = 0.5;
C = 0.5;
dx = Length/(Points);
dt = (C*dx)/U;

Matrix_Lax = zeros(t/dt, Points);
Matrix_Wen = zeros(t/dt, Points);

for i = 1:Points   %Apply initial solution for all points
    
    Matrix_Lax(1,i) = sin(2*pi*(i/Points));
    Matrix_Wen(1,i) = sin(2*pi*(i/Points));
    
end

%LAX Fredrichs method

for time = 2:(t/dt)
      
    Matrix_Lax(time,1) = -sin(2*pi*(dt*time));   %Left boundary condition
    
    for i =2:Points -1  %Lax fredrisch method
    
        Matrix_Lax(time, i) = 0.5*(Matrix_Lax(time - 1,i+1) + Matrix_Lax(time -1, i-1)) - (C/2)*(Matrix_Lax(time - 1, i+1) - Matrix_Lax(time - 1, i-1));

    end
    %For last point in the Vect, 2 point backwards is used

   Matrix_Lax(time, Points) = Matrix_Lax(time - 1, Points) - (C/2)*(3*Matrix_Lax(time - 1, Points) - 4*Matrix_Lax(time - 1, Points - 1) + Matrix_Lax(time -1, i -2));
   
end



%LAX Wendroff method

for time = 2:(t/dt)
      
    Matrix_Wen(time,1) = -sin(2*pi*(dt*time));   %Left boundary condition
    
    for i =2:Points -1  %Lax fredrisch method
    
        Matrix_Wen(time, i) = Matrix_Wen(time - 1, i) - (C/2)*(Matrix_Wen(time-1, i +1) - Matrix_Wen(time -1, i-1)) + ((C^2)/2)*(Matrix_Wen(time - 1, i- 1) - 2*Matrix_Wen(time - 1, i) + Matrix_Wen(time - 1, i +1));
        
    end
    %For last point in the Vect, 2 point backwards is used

   Matrix_Wen(time, Points) = Matrix_Wen(time - 1, Points) - (C/2)*(3*Matrix_Wen(time - 1, Points) - 4*Matrix_Wen(time - 1, Points - 1) + Matrix_Wen(time -1, i -2));
   
end

sin_vec = zeros(1, Points);
x = 0;
for x_vec = 1:Points
   
    sin_vec(x_vec) = sin(2*pi*(x-U*t));
     x = x + dx;
end

Error_Lax = abs(sin_vec - Matrix_Lax(end,:));
Error_Wen = abs(sin_vec - Matrix_Wen(end,:));

Norm_Lax = vecnorm(Error_Lax,1);
Norm_Wen = vecnorm(Error_Wen,1);

%Wendroof is better due to it including artificial diffusion


plot(Matrix_Lax(end,:),'b')
hold ON;
plot(Matrix_Wen(end,:),'r')
plot(sin_vec,'g')
% Create ylabel
ylabel({'Y'});
% Create xlabel
xlabel({'Points'});
legend('Lax Fredrich', 'Lax Wendroof', 'Analytical');
