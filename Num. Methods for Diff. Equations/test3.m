% Test project 3

%% Task 1.1
clc;
clear all;

% parameter values
N = 7; 
M = 110; 
tend =1; 
dx = 1/(N+1); 
dt = tend/M; 
xx = linspace(0, 1, N+2); 
xxshort = linspace(dx, 1-dx, N);        

tt = linspace(0, tend, M+1);            
[T,X]=meshgrid(tt,xx);
if N==1                                              
Tdx = toeplitz(-2)./(dx.^2);
else
Tdx = toeplitz([-2 1 zeros(1, N-2)])./(dx.^2);                       
end

% boundary conditions
g = @(x) -x.*(x-tend);
uold = g(xxshort);  
uold = uold(:);
unew = uold;
k=2;

% using eulerstep
for i = dt:dt:tend   
unew(:,k) = eulerstep(Tdx, uold, dt);       
uold = unew(:,k);
k = k+1;
end

% boundary conditions
unew = [zeros([1 M+1]); unew];
unew = [unew; zeros([1 M+1])];

%results
surf(X, T, abs(unew)); % varför detta absbelopp
dt/dx^2



%% Task 1.2

% Crank Nicolson method. trapezoidal rule

clc;
clear all;

% parameter values
N = 30; M = 100; tend =1; dx = 1/(N+1); dt = tend/M; xx = linspace(0, 1, N+2); xxshort = linspace(dx, 1-dx, N);  
mu = dt/(dx^2);
tt = linspace(0, tend, M+1);            
[T,X]=meshgrid(tt,xx);
if N==1                                              
Tdx = toeplitz(-2)./(dx.^2);
else
Tdx = toeplitz([-2 1 zeros(1, N-2)])./(dx.^2);                       
end

% boundary conditions
g = @(x) -x.*(x-tend);
uold = g(xxshort);  
uold = uold(:);
unew = uold;
k=2;

% using eulerstep
for i = dt:dt:tend   
unew(:,k) = TRstep(Tdx, uold, dt);       
uold = unew(:,k);
k = k+1;
end

% boundary conditions
unew = [zeros([1 M+1]); unew];
unew = [unew; zeros([1 M+1])];

%results
surf(X, T, unew); 
courant = dt/dx^2

%% Task 2.1 The advection equation

clc;
clear all;
N = 20; 
M = 100; 
tend =5;
dt = tend/M;
dx = 1/N;
g = @(x) exp(-100.*(x-0.5).^2);
xx = linspace(0, 1, N);
tt = linspace(0, tend, M+1); 
[T,X]=meshgrid(tt,xx);
u = g(xx);
u = u(:);
u1 = u;
a=1;
amu = a.*dt./dx;
k=1;
for i = dt:dt:tend
unew(:,k) = LaxWen(u, amu);
u = unew(:,k);
k = k+1;
end
unew = [u1, unew];
surf(T, X, unew); 
cfl = a*dt/dx
zlabel('u(t, x)')
xlabel('t')
ylabel('x')
title('u(t, x) with CFL number 1')
hold on
%%
norm = rms(unew);
plot(tt, norm);

%% Task 3 
clc;
clear all;
N = 100; 
M = 600; 
tend =1;
dt = tend/M;
dx = 1/N;
g = @(x) exp(-100.*(x-0.5).^2);
xx = linspace(0, 1, N); 
tt = linspace(0, tend, M+1); 
[T,X]=meshgrid(tt,xx);
u = g(xx);
u = u(:);
u1 = u;
a=70;
d = 0.1;

k=1;

for i = dt:dt:tend
unew(:,k) = convdiff(u, a, d, dt, dx);
u = unew(:,k);
k = k+1;
end

unew = [u1, unew];
mesh(T, X, unew); 

Pe = abs(a/d)
title('u(t, x)')
xlabel('t')
ylabel('x')
hold on
condition = Pe*dx

%% Task 4

clc;
clear all;
N = 300; 
M = 1000; 
tend =1;
dt = tend/M;
dx = 1/N;
g = @(x) exp(-100.*(x-0.5).^2);
xx = linspace(0, 1, N); 
tt = linspace(0, tend, M+1); 
[T,X]=meshgrid(tt,xx);
u = g(xx);
u = u(:);
u1 = u;
d = 0.01;

% create Tdx
sub = 1*ones(1, length(u)-1);
main = -2.*ones(1,length(u));
Tdx = diag(sub,-1) + diag(main) + diag(sub,1);
Tdx(1, length(u)) = 1;
Tdx(length(u), 1) = 1;
Tdx = Tdx./(dx^2);

% create A
A = eye(length(u))-d.*dt./2.*Tdx;

% solving loop
k=1;
for i = dt:dt:tend
unew(:,k) = A \(LW(u, dx, dt) + d*dt/2.*TRstep(Tdx, u, dt));
u = unew(:,k);
k = k+1;
end

unew = [u1, unew];
mesh(T, X, unew); 
title('u(t, x)')
xlabel('t')
ylabel('x')
hold on