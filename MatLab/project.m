%% FEM project first part
clc;

%% Stationary temperatur distribution
%%%% Parameters %%%%
coord = p';                            % coordinates from p
enod=t(1:3,:)';                        % nodes of elements
nelm=size(enod,1);                     % number of elements
nnod=size(coord,1);                    % number of nodes
dof=(1:nnod)';                         % degrees of freedom
ndof = length(dof);                    % number degrees of freedom
D = eye(2);                            % empty constitutive matrix
K = zeros(ndof);                       % empty stiffness matrix
fl = zeros(ndof,1);                    % empty load vector
fb = zeros(ndof,1);                    % empty load boundary vector
eq = zeros(1,nelm);                    % empty heat supply per unit vector
k = zeros(1,nelm);                     % empty thermal conductivity vector
Tinf = 18;                             % temperature far away
alphac = 40/(1000^2);                  % thermal transmittance
A = 10;                                % thickness
Q = (5 * 10^7)*10^(-9);                 % heat supply (Si)
%Q = 0.75*Q;                           % Netflix decreases video quality with 25%

%%%% Edof from mesh %%%%
dof_S=[(1:nnod)',(nnod+1:2*nnod)'];
for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

%%%% Convection segments %%%%
er = e([1 2 5],:);                     % Reduced e
conv_segments = [6 25 26 27];          % Chosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];    % edges with convection
    end
end

%%%% Coordinates and structure plot %%%%
[ex,ey]=coordxtr(edof,coord,dof,3);
%eldraw2(ex,ey,[1 4 0],edof(:,1)) % Plot of structure

%%%% Heat supply per unit vector - eq %%%%
for i = 1:length(t)
    if t(4,i) == 4
        eq(1,i) = Q;
    end
end

%%%% Thermal conductivity vector - k %%%%
for i = 1:length(t)
    if t(4,i) == 3 || t(4,i) == 1
        k(1,i) = 385*0.001;
    elseif t(4,i) == 2
        k(1,i) = 5*0.001;
    else
        k(1,i) = 149*0.001;
    end
end

%%%% Stiffnes matrix K and load vector fl %%%%
for elnr=1:nelm
    [Ke, fle] = flw2te(ex(elnr,:), ey(elnr,:), A, k(1,elnr)*D, eq(1,elnr));
    indx = edof(elnr,2:end);
    K(indx,indx) = K(indx,indx)+Ke;
    fl(indx) = fl(indx) + fle;
end

%%%% Boundary conditions bc, only used before convection is added %%%%
bc = [edges_conv(1,:) edges_conv(2,:)];
bc = unique(bc)';
bc(:,2) = 18;

%%%% Total load vector f %%%%
fb = A*assemfb(zeros(ndof,1),edges_conv,coord,alphac,Tinf);     % boundary load vector
f = fl + fb;

%%%% Total stiffness matrix %%%%
Kc = A*assemconv(zeros(ndof),edges_conv,coord,alphac);          % stiffness matrix from convection
K = K + Kc;


%%%% Solver A
%%%% Solve and extract %%%%
a = solveq(K, f);
et = extract(edof, a);
maxa = max(a); % max temperature 

%%%% Expand coordinates %%%%
exexp = [-flip(ex) ; ex];
eyexp = [flip(ey) ; ey];
etexp = [flip(et) ; et];

%%%% Plot result %%%%
patch(exexp',eyexp',etexp','EdgeColor','none')
title('Temperature distribution [C]')
colormap(autumn);
colorbar;
xlabel('x-position [mm]')
ylabel('y-position [mm]')
axis equal


%% Transient temperature distribution

%%%% Parameters %%%%
rho = zeros(1,nelm);                        % empty density vector
cp = zeros(1,nelm);                         % empty specific heat vector
C = zeros(ndof);                            % empty element matrix C
time = 20*60;                               % time from start in seconds
aold = 30.*ones(ndof,1);                    % start temperature at all nodes
deltat = 1;                                 % time steps
timevec = linspace(deltat, time, time)';    % time vector
anew = zeros(ndof,1);                       % temperature vector
max_a = zeros(time,1);                      % max temperature vector

%%%% Density rho and specific heat vector cp %%%%
for i = 1:length(t)
    if t(4,i) == 3 ||t(4,i)== 1
        rho(1,i) = 8930/10^9;
        cp(1,i) = 386;
    elseif t(4,i) == 2
        rho(1,i) = 2500/10^9;
        cp(1,i) = 1000;
    else
        rho(1,i) = 2530/10^9;
        cp(1,i) = 703;
    end
end

%%%% Element matrix - C %%%%
C = assemc(edof, ex, ey, C, cp, rho, A);

%%%% Implicit time stepping %%%%
for i = deltat:deltat:time
    anew = (C+deltat.*K)\(C*aold + deltat.*f);
    max_a(i,1) = max(anew);
    aold = anew;
end

a_max = max(max_a)      % max temperature


%%%% Extract %%%%
et = extract(edof, anew);

%%%% Expand coordinates %%%%
exexp = [-flip(ex) ; ex];
eyexp = [flip(ey) ; ey];
etexp = [flip(et) ; et];

%%%% Plot result %%%%
plot(timevec, max_a, '-');                          % plot max temp over time
%patch(exexp',eyexp',etexp','EdgeColor','none');    % plot tempdistribution
title('Maximum temperature, full quality ')
colormap(autumn);
xlabel('time [s]')
ylabel('temperature [C]')



%% Thermal expansion
clc;

%%%% Parameters %%%%
ndof2 = 2*ndof;                 % number degrees of freedom
Depsilon0 = zeros(3,nelm);      % epsilon0*D
E = [7 165 128].*10^9;          % Youngs modulus [Ag epoxy, Si, Cu]
v = [0.3 0.22 0.36];            % Poissons ratio [Ag epoxy, Si, Cu] *SKA DENNA VA MINUS?*
alpha = [40 2.6 17.6].*10^(-6); % Expansion coefficient [Ag epoxy, Si, Cu]
T0 = 30;                        % start temperature
deltaT = zeros(nelm, 1);        % vector for the elements difference in temperature
K2 = zeros(ndof2);              % stiffness matrix
D2 = zeros(3);                  % D matrix
ep = [2 A];                     % [p-type thickness]
f2 = zeros(ndof2, 1);           % thermal force vector
DAgEp = E(1)./((1+v(1))*(1-2*v(1))).*[1-v(1) v(1) 0; v(1) 1-v(1) 0 ; 0 0 0.5*(1-2*v(1))];   % D matrix for Ag Epoxy
DSi = E(2)./((1+v(2))*(1-2*v(2))).*[1-v(2) v(2) 0; v(2) 1-v(2) 0 ; 0 0 0.5*(1-2*v(2))];     % D matrix for Si
DCu = E(3)./((1+v(3))*(1-2*v(3))).*[1-v(3) v(3) 0; v(3) 1-v(3) 0 ; 0 0 0.5*(1-2*v(3))];     % D matrix for Cu
Seff_nod = zeros(nnod,1);       % Von Mises effective stress in the nodes
Seff_el = zeros(nelm,1);        % Von Mises effective stress in an element
sigma = zeros(nelm,3);          % stress

%%%% Delta T %%%%
for elnr = 1:nelm
    deltaT(elnr,1) = (a(t(1,elnr))+a(t(2,elnr))+a(t(3,elnr)))/3 - T0;
end

%%%% D2*epsilon0 %%%%
for elnr = 1:nelm
    if t(4,elnr) == 3 || t(4,elnr) == 1
        Depsilon0(1:2,elnr) =  alpha(3).*E(3).*deltaT(elnr)./(1-2.*v(3));    % Cu
    elseif t(4,elnr) == 2
        Depsilon0(1:2,elnr) =  alpha(1).*E(1).*deltaT(elnr)./(1-2.*v(1));   %Ag epoxy
    else
        Depsilon0(1:2,elnr) = alpha(2).*E(2).*deltaT(elnr)./(1-2.*v(2));    % Si
    end
end

%%%% Boundary conditions %%%%
% PCB board boundary
ed_segments = [12]; % Chosen boundary segments
edges_ed = [];
for i = 1:size(er,2)
    if ismember(er(3,i),ed_segments)
        edges_ed = [edges_ed er(1:2,i)];
    end
end

uniquefreedom1 = unique(edges_ed); % unique degrees of freedom
len = size(uniquefreedom1,1);
for i = 1:len
    for col = 2:7
        for row = 1:nelm
            if edof_S(row, col) == uniquefreedom1(i)
                if mod(col,2) == 0
                    uniquefreedom1 = [uniquefreedom1; edof_S(row,col+1)];   % fixing degrees of freedom
                elseif mod(col,2) ~= 0
                    uniquefreedom1 = [uniquefreedom1; edof_S(row,col-1)];   % fixing degrees of freedom
                    
                end
            end
        end
    end
end

% Left boundary in figure
ed_segments = [19 20 21 22]; % Chosen boundary segments
edges_ed = [];
for i = 1:size(er,2)
    if ismember(er(3,i),ed_segments)
        edges_ed = [edges_ed er(1:2,i)];
    end
end

uniquefreedom2 = unique(edges_ed);
len = size(uniquefreedom2,1);

for i = 1:len
    for col = 2:7
        for row = 1:nelm
            if edof_S(row, col) == uniquefreedom2(i)
                if mod(col,2) ~= 0
                    uniquefreedom2(i) = [uniquefreedom2; edof_S(row,col-1)];     % fixing degrees of freedom
                end
            end
        end
    end
end

% bc
uniquefreedom1 = unique(uniquefreedom1);
uniquefreedom2 = unique(uniquefreedom2);
uniquefreedomtot = [uniquefreedom1; uniquefreedom2];
bc = zeros(size(uniquefreedomtot,1), 2);

for i = 1:size(uniquefreedomtot,1)
    bc(i,1) = uniquefreedomtot(i);  % all fixated degrees of freedom
end

%%%% Stiffness matrix K2 and thermal force f2 %%%%
for elnr = 1: nelm
    if t(4,elnr) == 3 || t(4, elnr) == 1
        D2 = DCu;
    elseif t(4,elnr) == 2
        D2 = DAgEp;
    else
        D2 = DSi;
    end
    K2e = plante(ex(elnr,:),ey(elnr,:),ep,D2);
    f2e = plantf(ex(elnr,:),ey(elnr,:),ep,Depsilon0(:, elnr)');
    f2 = insert(edof_S(elnr, :), f2, f2e);
    indx = edof_S(elnr,2:end);
    K2(indx,indx) = K2(indx,indx)+K2e;
end

%%%% Solver %%%%
disp = solveq(K2, f2, bc);
ed = extract(edof_S, disp);

%%%% Expand coordinates %%%%
exexp = [-flip(ex) ; ex];
eyexp = [flip(ey) ; ey];

%%%% Calculate displaced coordinates %%%%
mag = 100; % Magnification (due to small deformations)
exd = ex + mag.*ed(:,1:2:end);
exdm =[-flip(exd) ; exd];
eyd = ey + mag.*ed(:,2:2:end);
eydm = [flip(eyd) ; eyd];
figure()
patch(exexp',eyexp',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
patch(exdm',eydm',[0 0 0],'FaceAlpha',0.3)
axis equal
title('Displacement field, full quality [Magnitude enhancement 100]')
xlabel('x-position [mm]')
ylabel('y-position [mm]')

%%
%%%% Von Mises %%%%
clc;

%%%% Stress sigma %%%%
for i = 1:nelm
    if t(4,i) == 3 || t(4, i) == 1
        D2 = DCu;
    elseif t(4,i) == 2
        D2 = DAgEp;
    else
        D2 = DSi;
    end
    es = plants(ex(i,:),ey(i,:),ep,D2,ed(i,:));
    sigma(i,:) = es;
end

sigma= sigma - Depsilon0';

%%%% Effective Von Mises stress %%%%
for i = 1:nelm
    if t(4,i) == 3 || t(4, i) == 1
        szz=v(3).*(sigma(i,1) + sigma(i,2)) - alpha(3)*E(3)*deltaT(i);
    elseif t(4,i) == 2
        szz=v(1).*(sigma(i,1) + sigma(i,2)) - alpha(1)*E(1)*deltaT(i);
    else
        szz=v(2).*(sigma(i,1) + sigma(i,2)) - alpha(2)*E(2)*deltaT(i);
    end
    Seff_el(i,1) = sqrt(sigma(i,1)^2 + sigma(i,2)^2 + szz^2 - sigma(i,1)*sigma(i,2) -sigma(i,1)*szz - sigma(i,2)*szz + 3*sigma(i,3)^2);
end

%%%% Maximum stress and coordinates %%%%
maxstress = max(Seff_el)
exmax = ex(elementnr,1)
eymax = ey(elementnr,1)

%%%% Node stress %%%% 
for i=1:2*nelm
    [c0,c1]=find(edof_S(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))./size(c0,1);
end

Seff_nod = extract(edof,Seff_nod);

%%%% Expand coordinates %%%%
exexp = [-flip(ex) ; ex];
eyexp = [flip(ey) ; ey];
Seff_nod = [flip(Seff_nod) ; Seff_nod];


%%%% Plot result %%%%
patch(exexp',eyexp',Seff_nod','EdgeColor','none');
title('Von Mises stress')
colormap(jet(5000));
colorbar;
xlabel('x-position [mm]')
ylabel('y-position [mm]')
axis equal
