%% 3 STEADY STATE TEMPERATURE USING 3 Central Diff, when 0.5*dx = dy
close all, clear all, clc
format long;

%Given data:

L = 0.4;     %[m] x-direction             
H = 0.2;     %[m] y-direction

prompt =  {'Input values for step size dx, dy = 0.5*dx'};
dlgtitle = 'Input';
definput = {'0.0125'};
dims = [1 50];
%Step size
Input_matrix = inputdlg(prompt, dlgtitle, dims, definput);
dx = (str2double(Input_matrix(1,1)));
dy = 0.5*dx;

%We got three different points with different x values and constant y value
%That we want to find temperature values for.

x1 = (0.1/dx)+1;
x2 = (0.2/dx)+1;
x3 = round((0.3/dx)+1);
y = (0.1/dy) +1;

            %[vertical(y), horizontal(x)]
            
%Amounts of points in the Rows/collums
M = (H/dy)+1; 
N = (L/dx)+1;
            
%   A)
%Makes the matrices for both step sizes
Matrix_1 = zeros(M, N);  

%Inputs the boundary condition (Corner pieces not correct = 0)
Top_BC = 20;
Bottom_BC = 10;
Left_BC = 30;
Right_BC = 100;
%BOUNDARY CONDITION APPLIED TO THE MATRIX
Matrix_1(M,2:N - 1) = Bottom_BC;
Matrix_1(1,2:N - 1) = Top_BC;
Matrix_1(2:M -1,1) = Left_BC;
Matrix_1(2:M -1,N) = Right_BC;

%Gives average values for corner pieces, just for the 2d contur plot
Matrix_1(1,1) = (Left_BC + Top_BC)/2;
Matrix_1(1,N) = (Right_BC + Top_BC)/2;
Matrix_1(M,1) = (Left_BC + Bottom_BC)/2;
Matrix_1(M,N) = (Right_BC + Bottom_BC)/2;

% Creates all the vectors and matrixes used to solve T_vector
Vect_1_b = zeros(1, ((M - 2)*(N - 2)));          

Matrix_1_A = zeros((M - 2)*(N - 2));

T_Vect_1 = ones(1, (M - 2)*(N - 2));


%3 Point difference with dx = dyc

%Corner pieces
for i = 1:1:length(Vect_1_b)


    if i == 1 %Left top corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i+1) = 1;
        Matrix_1_A(i, i+(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(1,2) + Matrix_1(2,1);  %Times four due to dy = d*0.5
    end
    
    
    if i == (N - 2)    %right top corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i-1) = 1;
        Matrix_1_A(i, i+(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(1,N -1) + Matrix_1(2,N);  %Times four due to dy = d*0.5
    end
    
    if i == ((N-2)*(M-3)) + 1    %left down corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i+1) = 1;
        Matrix_1_A(i, i-(N-2)) = 4;
        
        Vect_1_b(i) = Matrix_1(M -1 ,1) + 4*Matrix_1(M,2);  %Times four due to dy = d*0.5
    end
    
    if i == (N-2)*(M-2)      %Right down corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i-1) = 1;
        Matrix_1_A(i, i-(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(M,N -1) + Matrix_1(M - 1 ,N);  %Times four due to dy = d*0.5
    end
    

end


%Left most row of unknowns
k = 2;
for i =  (N-1):(N -2):(N-2)*(M-3)

    k = k +1;
    
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    Matrix_1_A(i,i-(N-2)) = 4;
    
    Vect_1_b(i) = Matrix_1(k,1);
    
end

%Top row of unkowns
k= 2;
for i = 2:1:(N-3)
    
    k = k +1;
    
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    
    Vect_1_b(i) = 4*Matrix_1(1,k);  %Times four due to dy = d*0.5
    
end
    
%central points of unkowns

for f = 0:1:(M - 5)
    
  
    for i = f*(N - 2)+(N):1:(f*(N - 2)+(N)) + (N - 5)
    
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    Matrix_1_A(i,i-(N-2)) = 4;
    
    
    end
end

%Right most row of unknowns
k = 2;
for i = 2*(N-2):(N-2):(N-2)*(M-3)
    k = k + 1;
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    Matrix_1_A(i,i-(N-2)) = 4;
    
    Vect_1_b(i) = Matrix_1(k,N);
    
end

%bottom row of unkowns
k = 2;
for i = ((N-2)*(M - 3))+2:1:((N-2)*(M-2))-1        
        k = k + 1;
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i-(N-2)) = 4;
    
     Vect_1_b(i) = 4*Matrix_1(M,k);  %Times four due to dy = d*0.5
end


T_vect_1 = Matrix_1_A \ -Vect_1_b';  %solve the equ system for Temp values

k = 0;
for b = 2:1:(M - 1)
    
    for i = 2:1:(N - 1)
        k = k+1;
    Matrix_1(b, i) =  T_vect_1(k);
    
    end
end %Fills the matrix with Temp values

%Steady state Temps at the sought of points x1, x2, x3 all at y heigth

X1 = ['At first point ', num2str(Matrix_1(y,x1)), ' Celcius'];
X2 = ['At second point ', num2str(Matrix_1(y,x2)), ' Celcius'];
X3 = ['At third point ', num2str(Matrix_1(y,x3)), ' Celcius'];

disp(X1)
disp(X2)
disp(X3)

contourf(Matrix_1)

%% Question 3: 5 Point Central differention Temp Steady state 0.5dx = dy
%When running this section a program that will ask for the dx will come up
%dy will be half of dx. It will give out the sought of temps in the command
%window and a figure with the contour of the matrix will be done.


close all, clear all, clc
format long;

%Given data:

L = 0.4;     %[m] x-direction             
H = 0.2;     %[m] y-direction

prompt =  {'Input values for step size dx, dy = 0.5*dx'};
dlgtitle = 'Input';
definput = {'0.0125'};
dims = [1 50];
%Step size
Input_matrix = inputdlg(prompt, dlgtitle, dims, definput);
dx = (str2double(Input_matrix(1,1)));
dy = 0.5*dx;      
 

%We got three different points with different x values and constant y value
%That we want to find temperature values for.

x1 = (0.1/dx)+1;
x2 = (0.2/dx)+1;
x3 = round((0.3/dx)+1);
y = (0.1/dy) +1;

            %[vertical(y), horizontal(x)]
            
%Amounts of points in the Rows/collums
M = (H/dy)+1; 
N = (L/dx)+1;
            
%   A)
%Makes the matrices for both step sizes
Matrix_1 = zeros(M, N);  

%Inputs the boundary condition (Corner pieces not correct = 0)
Top_BC = 20;
Bottom_BC = 10;
Left_BC = 30;
Right_BC = 100;
%BOUNDARY CONDITION APPLIED TO THE MATRIX
Matrix_1(M,2:N - 1) = Bottom_BC;
Matrix_1(1,2:N - 1) = Top_BC;
Matrix_1(2:M -1,1) = Left_BC;
Matrix_1(2:M -1,N) = Right_BC;

%Gives average values for corner pieces, just for the 2d contur plot
Matrix_1(1,1) = (Left_BC + Top_BC)/2;
Matrix_1(1,N) = (Right_BC + Top_BC)/2;
Matrix_1(M,1) = (Left_BC + Bottom_BC)/2;
Matrix_1(M,N) = (Right_BC + Bottom_BC)/2;

% Creates all the vectors and matrixes used to solve T_vector
Vect_1_b = zeros(1, ((M - 2)*(N - 2)));          

Matrix_1_A = zeros((M - 2)*(N - 2));

%3 Point difference with 0.5*dx = dy

%Corner pieces
for i = 1:1:length(Vect_1_b)


    if i == 1 %Left top corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i+1) = 1;
        Matrix_1_A(i, i+(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(1,2) + Matrix_1(2,1);  %Times four due to dy = d*0.5
    end
    
    
    if i == (N - 2)    %right top corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i-1) = 1;
        Matrix_1_A(i, i+(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(1,N -1) + Matrix_1(2,N);  %Times four due to dy = d*0.5
    end
    
    if i == ((N-2)*(M-3)) + 1    %left down corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i+1) = 1;
        Matrix_1_A(i, i-(N-2)) = 4;
        
        Vect_1_b(i) = Matrix_1(M -1 ,1) + 4*Matrix_1(M,2);  %Times four due to dy = d*0.5
    end
    
    if i == (N-2)*(M-2)      %Right down corner
        
        Matrix_1_A(i,i) = -10;
        Matrix_1_A(i,i-1) = 1;
        Matrix_1_A(i, i-(N-2)) = 4;
        
        Vect_1_b(i) = 4*Matrix_1(M,N -1) + Matrix_1(M - 1 ,N);  %Times four due to dy = d*0.5
    end
    

end

%Left most row of unknowns where 3 points central diff has to be used
k = 2;
for i =  (N-1):(N -2):(N-2)*(M-3)

    k = k +1;
    
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    Matrix_1_A(i,i-(N-2)) = 4;
    
    Vect_1_b(i) = Matrix_1(k,1);
    
end

%Top row of unkowns where 3 points central diff has to be used
k= 2;
for i = 2:1:(N-3)
    
    k = k +1;
    
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    
    Vect_1_b(i) = 4*Matrix_1(1,k);  %Times four due to dy = d*0.5
    
end
    
%central points of unkowns Where five point central can be used(Middle
%middle)

for f = 0:1:(M - 7)
    
  
    for i = f*(N-2)+((N*2)-1):1:f*(N-2)+((N*2)-1)+(N-7)
    
    Matrix_1_A(i,i) = -150;
    %Closes points
    Matrix_1_A(i,i+1) = 16;
    Matrix_1_A(i,i-1) = 16;
    Matrix_1_A(i,i+(N-2)) = 64;
    Matrix_1_A(i,i-(N-2)) = 64;
    %Second out points
    Matrix_1_A(i,i+2) = -1;
    Matrix_1_A(i,i-2) = -1;
    Matrix_1_A(i,i-(2*(N-2))) = -4;
    Matrix_1_A(i,i+(2*(N-2))) = -4;
    
    end
end

%Right most row of unknowns where 3 points central diff has to be used
k = 2;
for i = 2*(N-2):(N-2):(N-2)*(M-3)
    k = k + 1;
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+(N-2)) = 4;
    Matrix_1_A(i,i-(N-2)) = 4;
    
    Vect_1_b(i) = Matrix_1(k,N);
    
end

%bottom row of unkowns where 3 points central diff has to be used
k = 2;
for i = ((N-2)*(M - 3))+2:1:((N-2)*(M-2))-1        
        k = k + 1;
    Matrix_1_A(i,i) = -10;
    Matrix_1_A(i,i-1) = 1;
    Matrix_1_A(i,i+1) = 1;
    Matrix_1_A(i,i-(N-2)) = 4;
    
     Vect_1_b(i) = 4*Matrix_1(M,k);  %Times four due to dy = d*0.5
end


%Third row and collum have BC for the 5 points diff at the outer
%most points

%Corner pieces for the inner matrix where 5 point can be used but with BC
for i = 1:1:length(Vect_1_b)


    if i == N %Left top corner
        
        Matrix_1_A(i,i) = -150;
        %Closes points
        Matrix_1_A(i,i+1) = 16;
        Matrix_1_A(i,i-1) = 16;
        Matrix_1_A(i, i+(N-2)) = 64;
        Matrix_1_A(i, i-(N-2)) = 64;
        %Outer points
        Matrix_1_A(i,i+2) = -1;
        Matrix_1_A(i,i+((N-2)*2)) = -4;
     
        
        Vect_1_b(i) = -4*(Matrix_1(1,3) - Matrix_1(3,1));  %Times four due to dy = d*0.5
        
        
    end
    
    
    if i == N+(N-5)  %right top corner
        
        Matrix_1_A(i,i) = -150;
        %Closes pointa
        Matrix_1_A(i,i-1) = 16;
        Matrix_1_A(i,i+1) = 16;
        Matrix_1_A(i, i+(N-2)) = 64;
        Matrix_1_A(i, i-(N-2)) = 64;
        %Outer points
        Matrix_1_A(i,i-2) = -1;
        Matrix_1_A(i,i+((N-2)*2)) = -4;
        
        Vect_1_b(i) = -4*Matrix_1(1,N -2) - Matrix_1(3,N);  %Times four due to dy = d*0.5
    end
    
    
    if i == ((N-2)*(M-4)) + 2   %left down corner
        
        Matrix_1_A(i,i) = -150;
        %Closes pointa
        Matrix_1_A(i,i-1) = 16;
        Matrix_1_A(i,i+1) = 16;
        Matrix_1_A(i, i+(N-2)) = 64;
        Matrix_1_A(i, i-(N-2)) = 64;
        %Outer points
        Matrix_1_A(i,i+2) = -1;
        Matrix_1_A(i,i-((N-2)*2)) = -4;
        
        Vect_1_b(i) = -Matrix_1(M-2,3) - 4*Matrix_1(M,3);  %Times four due to dy = d*0.5
    end
    
    if i == (N-2)*(M-3) -1      %Right down corner
        
        Matrix_1_A(i,i) = -150;
        %Closes pointa
        Matrix_1_A(i,i-1) = 16;
        Matrix_1_A(i,i+1) = 16;
        Matrix_1_A(i, i+(N-2)) = 64;
        Matrix_1_A(i, i-(N-2)) = 64;
        %Outer points
        Matrix_1_A(i,i-2) = -1;
        Matrix_1_A(i,i-((N-2)*2)) = -4;
        
        Vect_1_b(i) = -Matrix_1(M-2,N) - 4*Matrix_1(M,N-2);  %Times four due to dy = d*0.5
    end
    

end


%Top row where 5 point can be used but have boundary condition in the outer
%points:

k = 3;
for i = N+1:1:N+1+(N-7)
    
    k = k + 1;
    Matrix_1_A(i,i) = -150;
    %Close points
    Matrix_1_A(i,i-1) = 16;
    Matrix_1_A(i,i+1) = 16;
    Matrix_1_A(i,i-(N-2)) = 64;
    Matrix_1_A(i,i+(N-2)) = 64;
    %Outer points
    Matrix_1_A(i,i+2) = -1;
    Matrix_1_A(i,i-2) = -1;
    Matrix_1_A(i,i+((N-2)*2)) = -4;
    
    Vect_1_b(i) = -4*Matrix_1(1,k);  %Times four due to dy = d*0.5

end

%Left side row where 5 point can be used but have boundary condition in the
%outer points

k = 3;
for i = ((2*N)-2):N-2:((2*N)-2)+(M-7)*(N-2)
    
    k = k + 1;
    Matrix_1_A(i,i) = -150;
    %Close points
    Matrix_1_A(i,i-1) = 16;
    Matrix_1_A(i,i+1) = 16;
    Matrix_1_A(i,i-(N-2)) = 64;
    Matrix_1_A(i,i+(N-2)) = 64;
    %Outer points
    Matrix_1_A(i,i+2) = -1;
    Matrix_1_A(i,i-((N-2)*2)) = -4;
    Matrix_1_A(i,i+((N-2)*2)) = -4;
    
    Vect_1_b(i) = -Matrix_1(k,1);

end

%right side row where 5 point can be used but have boundary condition in
%one point
k = 3;
for i = (((N-2)*3)-1):(N-2):((((N-2)*3)-1)+(M-7)*(N-2))
    k = k + 1;
    Matrix_1_A(i,i) = -150;
    %Close points
    Matrix_1_A(i,i-1) = 16;
    Matrix_1_A(i,i+1) = 16;
    Matrix_1_A(i,i-(N-2)) = 64;
    Matrix_1_A(i,i+(N-2)) = 64;
    %Outer points
    Matrix_1_A(i,i-2) = -1;
    Matrix_1_A(i,i-((N-2)*2)) = -4;
    Matrix_1_A(i,i+((N-2)*2)) = -4;
    
    Vect_1_b(i) = -Matrix_1(k,N);

end

%Bottom row where 5 point can be used but have boundary contitons in one
%point

k = 3;
for i = (((N-2)*(M-4)) + 3):1:(((N-2)*(M-4)) + 3)+(N-7)
    
    k = k + 1;
    Matrix_1_A(i,i) = -150;
    %Close points
    Matrix_1_A(i,i-1) = 16;
    Matrix_1_A(i,i+1) = 16;
    Matrix_1_A(i,i-(N-2)) = 64;
    Matrix_1_A(i,i+(N-2)) = 64;
    %Outer points
    Matrix_1_A(i,i+2) = -1;
    Matrix_1_A(i,i-2) = -1;
    Matrix_1_A(i,i-((N-2)*2)) = -4;
    
    Vect_1_b(i) = -4*Matrix_1(M,k);  %Times four due to dy = d*0.5
    
end

T_Vect_1 = Matrix_1_A \ -Vect_1_b';  %solve the equ system for Temp values

k = 0;
for b = 2:1:(M - 1)
    
    for i = 2:1:(N - 1)
        k = k+1;
    Matrix_1(b, i) =  T_Vect_1(k);
    
    end
end %Fills the matrix with Temp values

%Steady state Temps at the sought of points x1, x2, x3 all at y heigth

X1 = ['At first point ', num2str(Matrix_1(y,x1)), ' Celcius'];
X2 = ['At second point ', num2str(Matrix_1(y,x2)), ' Celcius'];
X3 = ['At third point ', num2str(Matrix_1(y,x3)), ' Celcius'];

disp(X1)
disp(X2)
disp(X3)

contourf(Matrix_1)




