% Main function
%% Task 1.1
clear all;
clc;

alpha = 0;
beta = 3;
L = 1;
N = 3;
dx = L/(N+1);
xvec = 0:dx:L;
f = @(x) 2+6*x;             % andraderivata till y=x+x^2+x^3
yexakt = @(x1) x1+x1^2+x1^3; 
for i = 1:N
  fvec(i) = f(xvec(i+1));   % diskretisering, använder i + 1 sådan att inte x0 kommer med
end
fvec = fvec(:);             % gör om till kolumnvektor
ynew = twopBVP(fvec, alpha, beta, L, N);

plot(xvec, ynew);           %plot för jämförelse mellan lösningarna
hold on;
fplot(yexakt, [0, L]);

%% Task 1.1 plotta error
clear all;
clc;

alpha = 5;
beta = 2;
L = 1;
f = @(x) 12*(beta-alpha)*x^2/(L^4);                                 % andraderivatan
yexakt = @(x1) (beta-alpha)*x1^4/(L^4) + alpha; 

for N = 1:100
    dx(N) = L/(N+1);                            % nytt dx för varje loop
    xvec = dx(N):dx(N):L-dx(N);                 % ny x-vektor för varje loop
    
    for i = 1:N
    fvec(i) = f(xvec(i));                       % diskretisering av f
    exakt(i) = yexakt(xvec(i));                 % diskretisering av exakta lösningen
    end
    
    fvec = fvec(:);
    exakt = exakt(:);
    exakt = [alpha; exakt; beta];
    ynew = twopBVP(fvec, alpha, beta, L, N);    % får ut hela ynew inklusive BV 
    globerr (N) = rms(ynew - exakt) ;           % rms (tar värde för hela vektorernas skillnad)
    exakt = zeros;                              % denna behövs då det sista elementet annars följer med.

end


loglog(dx, globerr); %varför får jag ut så "bra" värden?
hold on;
grid on;
title ('Error plot rms');
xlabel('\Deltax');
ylabel('error');
legend('error', '');

%% Task 1.2
% Solving Mbis = q(x)
clear all;clc;
alpha = 0;beta = 0;L = 10;N = 999;
qvec = -50*10^3;
dx = L/(N+1);
xvec = dx:dx:L-dx;                                 
for i = 1:N-1
qvec = [ qvec ; -50*10^3];                          
end
M = twopBVP(qvec, alpha, beta, L, N);
M(1) = [];
M(length(M)) = [];                                  
E = 1.9.*10.^(11);
I = @(x) 10.^(-3).*(3-2.*(cos(pi.*x./L)).^12);
for i = 1:N
   Ivec(i) = I(xvec(i));                           
end
Ivec = Ivec(:);
ubis = M./(Ivec.*E);                                 
u = 1000.*twopBVP(ubis, alpha, beta, L, N);         
xvec = [0 xvec L];                                  
plot(xvec, u);
hold on;
title ('Deflection u as a function of x');
ylabel('u(x) [mm]');
xlabel('x [m]');
plot(xvec(length(u)/2+0.5), u(length(u)/2+0.5), '*')
hold on
format long
u(length(u)/2+0.5)                                  
text(5,-11, 'u_{min}= -11.741059 mm')


%% Task 2.1
clear all;
clc;

alpha = 0;
beta = 0;
L = 1;
N = 499;

for i = 1:N
eigvalexakt(i) = -(i.*pi).^2;
end
eigvalexakt = eigvalexakt(:);
[eigenvec, eigenval] = sturm(L, N);
eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend'); % sortering sådan att den med minsta storleken kommer först


error = eigenval - eigvalexakt;
loglog([1 2 3], [error(1) error(2) error(3)]);
hold on

[~, eigenvalny] = sturm(L, 3);
eigenvalny = diag(eigenvalny);
eigenvalny = sort(eigenvalny, 'descend');
for i = 1:3
eigvalexaktny(i) = -(i.*pi).^2;
end
eigvalexaktny = eigvalexaktny(:);
errorny = eigenvalny - eigvalexaktny;
loglog([1 2 3], errorny, 'magenta');
hold on
grid on
title ('Error for each eigenvalue as a function of a specific N');
ylabel('error for eigenvalue(N)');
xlabel('N');

%% Task 2.1 

clear all;
clc;

N = 499;
L=1;

for i = 1:N
for j = 1:i
eigvalexakt(j) = -(j.*pi).^2;                   % gör en vektor av de exakta egenvärdena i varje loop
end 
eigvalexakt = eigvalexakt(:);                   % gör om den till en kolonnvektor för att senare matcha approximerade   

[~, eigenval] = sturm(L, i);
eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');               % sortering sådan att den med minsta storleken kommer först

error1(i) = rms(eigenval(1) - eigvalexakt(1));      %sparar rms-felet för första egenvärdet i en vektor
if i>1
error2(i) = rms(eigenval(2) - eigvalexakt(2));      %sparar rms-felet för andra egenvärdet i en vektor
end
if i>2
  error3(i) = rms(eigenval(3) - eigvalexakt(3));    %sparar rms-felet för tredje egenvärdet i en vektor
end
end

loglog(1:N, error1);
hold on;
loglog(1:N, error2, 'magenta');
hold on;
loglog(1:N, error3, 'green');
grid on;
title ('Error versus N');
ylabel('error(N)');
xlabel('N');
legend('First eigenvalue', 'Second eigenvalue', 'Third eigenvalue');

%% Task 2.1 eigenvalues
clear all;
clc;

N = 499;
L=1;

[~, eigenval] = sturm(L, N);
eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
display('First eigenval');
eigenval(1)
display('Second eigenval');
eigenval(2)
display('Third eigenval');
eigenval(3)

%% Task 2.1 eigenmodes
clc;
clear all;
N = 499;
L = 1;
dx = L./(N+1);
xvec = 0:dx:L;
[eigenvec, ~] = sturm(L, N);
eigenvec1 = -[0 ; eigenvec(:,N) ; 0];           % varför hamnar de upp och ner ? två av dem
eigenvec2 = -[0 ; eigenvec(:,N-1) ; 0];
eigenvec3 = [0 ; eigenvec(:,N-2) ; 0];

plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;
plot(xvec, eigenvec3, 'blue'); hold on;

legend('First eigenfunction', 'Second eigenfunction', 'Third eigenfunction');
title ('Eigenfunctions');
ylabel('y(x)');
xlabel('x');


%% Task 22 adderar Ek probabilities


clc;
clear all;
N = 499;
L = 1;
dx = L./(N+1);
xvec = 0:dx:L;
V=0;
[eigenvec, eigenval] = sturm(L, N, V);
eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (80.*([0 ; eigenvec(:,N) ; 0])).^2 + norm(eigenval(1));        
eigenvec2 = (80.*([0 ; eigenvec(:,N-1) ; 0])).^2 + norm(eigenval(2));                            
eigenvec3 = (80.*([0 ; eigenvec(:,N-2) ; 0])).^2 + norm(eigenval(3));
eigenvec4 = (80.*([0 ; eigenvec(:,N-3) ; 0])).^2 + norm(eigenval(4));        
eigenvec5 = (80.*([0 ; eigenvec(:,N-4) ; 0])).^2 + norm(eigenval(5));                            
eigenvec6 = (80.*([0 ; eigenvec(:,N-5) ; 0])).^2 + norm(eigenval(6));

plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;

legend('|\Psi´_1(x)|^2 + E_1', '|\Psi´_2(x)|^2 + E_2', '|\Psi´_3(x)|^2 + E_3', '|\Psi´_4(x)|^2 + E_4', '|\Psi´_5(x)|^2 + E_5', '|\Psi´_6(x)|^2 + E_6');
title ('Normalized probability density functions');
ylabel('|\Psi´_k(x)|^2 + E_k');
xlabel('x');

%% Task 22 adderar Ek vågfunktioner


clc;
clear all;
N = 499;
L = 1;
dx = L./(N+1);
xvec = 0:dx:L;
V=0;
[eigenvec, eigenval] = sturm(L, N, V);
eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-250.*([0 ; eigenvec(:,N) ; 0])) + norm(eigenval(1));        
eigenvec2 = (-250.*([0 ; eigenvec(:,N-1) ; 0])) + norm(eigenval(2));                            
eigenvec3 = (250.*([0 ; eigenvec(:,N-2) ; 0])) + norm(eigenval(3));
eigenvec4 = (250.*([0 ; eigenvec(:,N-3) ; 0])) + norm(eigenval(4));        
eigenvec5 = (-250.*([0 ; eigenvec(:,N-4) ; 0])) + norm(eigenval(5));                            
eigenvec6 = (-250.*([0 ; eigenvec(:,N-5) ; 0])) + norm(eigenval(6));

plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
legend('\Psi´_1(x) + E_1', '\Psi´_2(x) + E_2', '\Psi´_3(x) + E_3', '\Psi´_4(x) + E_4', '\Psi´_5(x) + E_5', '\Psi´_6(x) + E_6');
title ('Normalized wavefunctions');
ylabel('\Psi´(x)_k + E_k');
xlabel('x');

%% Task 2.2 First potential V wavefunctions
clear all;
clc;
V = @(x) 700.*(0.5 - abs(x - 0.5));
L=1;
N=499;
dx = L./(N+1);
xvec = dx:dx:L-dx;
Vmat = diag(V(xvec));
xvec = [0 xvec L];
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-350.*([0 ; eigenvec(:,N) ; 0])) + norm(eigenval(1)) ; 
eigenvec2 = (-350.*([0 ; eigenvec(:,N-1) ; 0])) + norm(eigenval(2));                            
eigenvec3 = (350.*([0 ; eigenvec(:,N-2) ; 0])) + norm(eigenval(3));
eigenvec4 = (350.*([0 ; eigenvec(:,N-3) ; 0])) + norm(eigenval(4));        
eigenvec5 = (-350.*([0 ; eigenvec(:,N-4) ; 0])) + norm(eigenval(5));                            
eigenvec6 = (-350.*([0 ; eigenvec(:,N-5) ; 0])) + norm(eigenval(6));
plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
legend('\Psi´_1(x) + E_1', '\Psi´_2(x) + E_2', '\Psi´_3(x) + E_3', '\Psi´_4(x) + E_4', '\Psi´_5(x) + E_5', '\Psi´_6(x) + E_6');
title ('Normalized wavefunctions');
ylabel('\Psi´(x)_k + E_k');
xlabel('x');

%% Task 2.2 First potential V probabilities
clear all;
clc;
V = @(x) 700.*(0.5 - abs(x - 0.5));
L=1;
N=499;
dx = L./(N+1);
xvec = dx:dx:L-dx;
Vmat = diag(V(xvec));
xvec = [0 xvec L];
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-100.*([0 ; eigenvec(:,N) ; 0])).^2 + norm(eigenval(1)) ; 
eigenvec2 = (-100.*([0 ; eigenvec(:,N-1) ; 0])).^2 + norm(eigenval(2));                            
eigenvec3 = (100.*([0 ; eigenvec(:,N-2) ; 0])).^2 + norm(eigenval(3));
eigenvec4 = (100.*([0 ; eigenvec(:,N-3) ; 0])).^2 + norm(eigenval(4));        
eigenvec5 = (100.*([0 ; eigenvec(:,N-4) ; 0])).^2 + norm(eigenval(5));                            
eigenvec6 = (100.*([0 ; eigenvec(:,N-5) ; 0])).^2 + norm(eigenval(6));
plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
legend('|\Psi´_1(x)|^2 + E_1', '|\Psi´_2(x)|^2 + E_2', '|\Psi´_3(x)|^2 + E_3', '|\Psi´_4(x)|^2 + E_4', '|\Psi´_5(x)|^2 + E_5', '|\Psi´_6(x)|^2 + E_6');
title ('Normalized probability density functions');
ylabel('|\Psi´_k(x)|^2 + E_k');
xlabel('x');

%% Task 2.2 second V wavefunctions
clear all;
clc;
V = @(x) 800.*(sin(pi.*x)).^2;
L=1;
N=499;
dx = L./(N+1);
xvec = dx:dx:L-dx;
Vmat = diag(V(xvec));
xvec = [0 xvec L];
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-350.*([0 ; eigenvec(:,N) ; 0])) + norm(eigenval(1)) ; 
eigenvec2 = (-350.*([0 ; eigenvec(:,N-1) ; 0])) + norm(eigenval(2));                            
eigenvec3 = (350.*([0 ; eigenvec(:,N-2) ; 0])) + norm(eigenval(3));
eigenvec4 = (350.*([0 ; eigenvec(:,N-3) ; 0])) + norm(eigenval(4));        
eigenvec5 = (-350.*([0 ; eigenvec(:,N-4) ; 0])) + norm(eigenval(5));                            
eigenvec6 = (-350.*([0 ; eigenvec(:,N-5) ; 0])) + norm(eigenval(6));
plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
legend('\Psi´_1(x) + E_1', '\Psi´_2(x) + E_2', '\Psi´_3(x) + E_3', '\Psi´_4(x) + E_4', '\Psi´_5(x) + E_5', '\Psi´_6(x) + E_6');
title ('Normalized wavefunctions');
ylabel('\Psi´(x)_k + E_k');
xlabel('x');

%% Task 2.2 second potential V probabilities
clear all;
clc;
V = @(x) 800.*(sin(pi.*x)).^2;
L=1;
N=499;
dx = L./(N+1);
xvec = dx:dx:L-dx;
Vmat = diag(V(xvec));
xvec = [0 xvec L];
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-100.*([0 ; eigenvec(:,N) ; 0])).^2 + norm(eigenval(1)) ; 
eigenvec2 = (-100.*([0 ; eigenvec(:,N-1) ; 0])).^2 + norm(eigenval(2));                            
eigenvec3 = (100.*([0 ; eigenvec(:,N-2) ; 0])).^2 + norm(eigenval(3));
eigenvec4 = (100.*([0 ; eigenvec(:,N-3) ; 0])).^2 + norm(eigenval(4));        
eigenvec5 = (100.*([0 ; eigenvec(:,N-4) ; 0])).^2 + norm(eigenval(5));                            
eigenvec6 = (100.*([0 ; eigenvec(:,N-5) ; 0])).^2 + norm(eigenval(6));
plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
legend('|\Psi´_1(x)|^2 + E_1', '|\Psi´_2(x)|^2 + E_2', '|\Psi´_3(x)|^2 + E_3', '|\Psi´_4(x)|^2 + E_4', '|\Psi´_5(x)|^2 + E_5', '|\Psi´_6(x)|^2 + E_6');
title ('Normalized probability density functions');
ylabel('|\Psi´_k(x)|^2 + E_k');
xlabel('x');

%% Task 2.2 Dirac potential V probabilities
clear all;clc;
L=1;N=499;V = zeros(1,N,'uint32');
V(round(N/3)) = 100000;V(2*round(N/3)) = 100000;
dx = L./(N+1);xvec = 0:dx:L;
Vmat = diag(double(V));
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);eigenval = sort(eigenval, 'descend');eigenval = eigenval(:);

eigenvec1 = (-100.*([0 ; eigenvec(:,N) ; 0])).^2 + norm(eigenval(1)) ; 
eigenvec2 = (-100.*([0 ; eigenvec(:,N-1) ; 0])).^2 + norm(eigenval(2));                            
eigenvec3 = (100.*([0 ; eigenvec(:,N-2) ; 0])).^2 + norm(eigenval(3));
eigenvec4 = (100.*([0 ; eigenvec(:,N-3) ; 0])).^2 + norm(eigenval(4));        
eigenvec5 = (100.*([0 ; eigenvec(:,N-4) ; 0])).^2 + norm(eigenval(5));                            
eigenvec6 = (100.*([0 ; eigenvec(:,N-5) ; 0])).^2 + norm(eigenval(6));
plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, eigenvec4); hold on;
plot(xvec, eigenvec5); hold on;                     
plot(xvec, eigenvec6); hold on;
plot(xvec, [0 double(V)./100 0 ], '-');hold on;
fplot(-eigenval(1), [0, 1], 'magenta');hold on;
fplot(-eigenval(2), [0, 1], 'green'); hold on;
fplot(-eigenval(3), [0, 1], 'blue'); hold on;
fplot(-eigenval(4), [0, 1]); hold on;
fplot(-eigenval(5), [0, 1]); hold on;
fplot(-eigenval(6), [0, 1]); hold on;

legend('|\Psi´_1(x)|^2 + E_1', '|\Psi´_2(x)|^2 + E_2', '|\Psi´_3(x)|^2 + E_3', '|\Psi´_4(x)|^2 + E_4', '|\Psi´_5(x)|^2 + E_5', '|\Psi´_6(x)|^2 + E_6', 'Potential V(x)');
title ('Normalized probability density functions');
ylabel('|\Psi´_k(x)|^2 + E_k');xlabel('x');

%% Task 2.2 Dirac potential V wave function
clear all;
clc;
L=1;
N=499;
V = zeros(1,N,'uint32');
V(round(N/3)) = 100000;
V(2*round(N/3)) = 100000;
dx = L./(N+1);
xvec = 0:dx:L;
Vmat = diag(double(V));
[eigenvec, eigenval] = sturm(L, N, Vmat);

eigenval = diag(eigenval);
eigenval = sort(eigenval, 'descend');
eigenval = eigenval(:);

eigenvec1 = (-100.*([0 ; eigenvec(:,N) ; 0])) + norm(eigenval(1)) ; 
eigenvec2 = (-100.*([0 ; eigenvec(:,N-1) ; 0])) + norm(eigenval(2));                            
eigenvec3 = (100.*([0 ; eigenvec(:,N-2) ; 0])) + norm(eigenval(3));

plot(xvec, eigenvec1, 'magenta'); hold on;
plot(xvec, eigenvec2, 'green'); hold on;                     
plot(xvec, eigenvec3, 'blue'); hold on;
plot(xvec, [0 double(V) 0 ], '-');hold on;
legend('\Psi´_1(x) + E_1', '\Psi´_2(x) + E_2', '\Psi´_3(x) + E_3', 'Potential V(x)');
title ('Normalized wave functions and potential V(x)');
ylabel('\Psi´_k(x) + E_k');
xlabel('x');
%%
clc;
fplot(@(x) 700.*(0.5 -abs(0.5-x)), [0, 1]);
