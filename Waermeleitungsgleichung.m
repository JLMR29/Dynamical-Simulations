clc
clear

% Definition des Meschgitters:
Mesh = ones(200);

% Anfangsbedingungen:
T_0 = 273.15;  %K
Mesh = T_0.*Mesh;

% Randbedingungen: konstante Temperatur an den Kanten
T_rand = 1000; %K
Mesh(1,:) = T_rand;
Mesh(200,:) = T_rand;
Mesh(:, 1) =T_rand;
Mesh(:, 200) = T_rand;
MeshNew = Mesh;
dudt = 0;
a=0.001;   % a= lambda/(rho*cp) --> if too high the calculation becomes unstable.
dx = 0.01;  %m
dy = 0.01;  %m
dt=0.01;    %s


for k=1:1000
for i=2:199
    for j=2:199
        dTdt = a.*((Mesh(i+1,j)-2.*Mesh(i,j)+Mesh(i-1,j))/(dx*dx) + (Mesh(i,j+1)-2.*Mesh(i,j)+Mesh(i,j-1))/(dy*dy));
        MeshNew(i,j) = Mesh(i,j)+dt.*dTdt;
    end
end
time = append(string(k*dt), " [s]");
Mesh = MeshNew;
imagesc(MeshNew);
title(time)
c = colorbar;
xlabel("x [cm]")
ylabel("y [cm]")
c.Label.String = "Temperature [K]";
pause(0.001);
end

