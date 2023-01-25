% This scripts solves the 2D heath transfer equation
% The method used for solving is backward euler (implicit)
clear
clc

% Generating the mesh:
NumOfpoints = 200;
Mesh = ones(NumOfpoints);

% Physical parameters:
a = 0.001;   % lambda/(rho*cp) 
x = 1;  % m
y = 1;  % m
dx = x/NumOfpoints;  %m
dy = y/NumOfpoints;  %m
T = 100; %s
dt = 0.001;
T_sim = T/dt;
display(T_sim)

% Boundary and initial conditions:
T_0 = 273.15;   %K
Mesh = T_0.*Mesh;
T_boundary_0 = 600;     %K
Mesh(1,:)=T_boundary_0;
Mesh(NumOfpoints,:)=T_boundary_0;
Mesh(:,1)=T_boundary_0;
Mesh(:,NumOfpoints)=T_boundary_0;
MeshNew = Mesh;

% For Thomas Algorithm:
A = ones(1,NumOfpoints-2);
B = ones(1,NumOfpoints-2);
C = ones(1,NumOfpoints-2);
D = ones(1,NumOfpoints-2);
TriDMatrix = zeros(NumOfpoints-2);
SolVector = zeros(1,NumOfpoints-2);
% Backwards difference using the ADI-scheme:
for t=1:T_sim
    % first plot the current mesh
    imagesc(Mesh)
    time = append("time= ",string((t-1).*dt));
    title(append(time," [s]"))
    xticklabels = 0:0.1:1;
    xticks = linspace(1, size(Mesh, 2), numel(xticklabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    yticklabels = 0:0.1:1;
    yticks = linspace(1, size(Mesh, 1), numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
    c = colorbar;
    xlabel("x [m]")
    ylabel("y [m]")
    c.Label.String = "Temperature [K]";
    pause(0.001)
    

    % Now actualize in x-Direction: 
    for j=2:(NumOfpoints-1)
        for i=2:(NumOfpoints-1)
            r = dt.*a/(dx.*dx);
            B = (1+2.*r);
            C = -r;
            A = -r;
            
            if i>2
                TriDMatrix(i-1,i-2)=A;
            end
            TriDMatrix(i-1,i-1) = B;

            if i<NumOfpoints-1
                TriDMatrix(i-1,i) = C;
            end

            if i==2
                D = (1-2.*r).*Mesh(i,j)+r.*(Mesh(i,j-1) + Mesh(i,j+1)) + r.*MeshNew(i-1,j);
            elseif i==(NumOfpoints-1)
                D = (1-2.*r).*Mesh(i,j)+r.*(Mesh(i,j-1) + Mesh(i,j+1)) + r.*MeshNew(i+1,j);
            else
                D = (1-2.*r).*Mesh(i,j)+r.*(Mesh(i,j-1)+Mesh(i,j+1));
            end
        SolVector(i-1) = D;
        
        end
        RightSide = transpose(SolVector);
        X = linsolve(TriDMatrix, RightSide);
        
        %Insert the solution in x-direction to the new mesh:
        MeshNew(j,2:(NumOfpoints-1)) = transpose(X);         % Boundary conditions (1 and NumOfpoints) are already implemented;
    end

    
    % Now in Y-Direction:
    % In this part of the code, Mesh takes the role of the new Mesh and
    % MeshNew is now the one getting actualized
    for i=2:(NumOfpoints-1)
        for j=2:(NumOfpoints-1)
            r = dt.*a/(dx.*dx);
            B = (1+2.*r);
            C = -r;
            A = -r;
            
            if j>2
                TriDMatrix(j-1,j-2)=A;
            end
            TriDMatrix(j-1,j-1) = B;

            if j<NumOfpoints-1
                TriDMatrix(j-1,j) = C;
            end

            if j==2
                D = (1-2.*r).*MeshNew(i,j)+r.*(MeshNew(i-1,j) + MeshNew(i+1,j)) + r.*Mesh(i,j-1);
            elseif j==(NumOfpoints-1)
                D = (1-2.*r).*MeshNew(i,j)+r.*(MeshNew(i-1,j) + MeshNew(i+1,j)) + r.*Mesh(i,j+1);
            else
                D = (1-2.*r).*MeshNew(i,j)+r.*(MeshNew(i-1,j) + MeshNew(i+1,j));
            end
        SolVector(j-1) = D;
        
        end
        RightSide = transpose(SolVector);
        Y = linsolve(TriDMatrix, RightSide);
        
        %Insert the solution in y-direction into the new mesh:
        Mesh(2:(NumOfpoints-1),i) = Y;         % Boundary conditions (1 and NumOfpoints) are already implemented;
    end

end


