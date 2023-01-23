clear
clc
% System constans:
m1 = 10;
m2 = 10;
m3 = 10;
k = 10;
deq = 1;
dt = 0.01;
T = 1000;

% starting conditions:
p1_prime = [0 0];
p2_prime = [0 0];
p3_prime = [0 0];
p1 = [0.2 0.5];
p2 = [0 1];
p3 = [1 0];

% vectors for saving data:
Position1 = zeros(T,2);
Position2 = zeros(T,2);
Position3 = zeros(T,2);
Time = zeros(T,1);
Distance12 = zeros(T,1);




for t=0:(T-1)
    r12=p2-p1;
    r23=p3-p2;
    r31=p1-p3;
    f12 = (-1)*k*(norm(r12)-deq)/norm(r12);
    f23 = (-1)*k*(norm(r23)-deq)/norm(r23);
    f31 = (-1)*k*(norm(r31)-deq)/norm(r31);
    F12 = f12.*r12;
    F23 = f23.*r23;
    F31 = f31.*r31;
    F21 = (-1).*F12;
    F32 = (-1).*F23;
    F13 = (-1).*F31;
    display(t)
   
    NetF1 = F21 + F31;
    NetF2 = F12 + F32;
    NetF3 = F23 + F13;
    p1_dprime = (1/m1).*NetF1;
    p2_dprime = (1/m2).*NetF2;
    p3_dprime = (1/m3).*NetF3;
    

    p1_pnplus1 = p1_prime + dt.*p1_dprime;
    p1_nplus1 = p1 + dt.*p1_prime+((dt*dt)/2).*p1_dprime;
    
    p2_pnplus1 = p2_prime + dt.*p2_dprime;
    p2_nplus1 = p2 + dt.*p2_prime+((dt*dt)/2).*p2_dprime;

    p3_pnplus1 = p3_prime + dt.*p3_dprime;
    p3_nplus1 = p3 + dt.*p3_prime+((dt*dt)/2).*p3_dprime;

    scatter(p1(1,1),p1(1,2),"filled","g"); hold on 
    scatter(p2(1,1),p2(1,2),"filled","r"); 
    scatter(p3(1,1),p3(1,2),"filled","b");
    


    pause(0.001);

    %saving data:
    Position1(t+1,1) = p1(1,1);
    Position1(t+1,2) = p1(1,2);
    Position2(t+1,1) = p2(1,1);
    Position2(t+1,2) = p2(1,2);
    Position3(t+1,1) = p3(1,1);
    Position3(t+1,2) = p3(1,2);


    %actualizing data:
    p1 = p1_nplus1;
    p2 = p2_nplus1;
    p3 = p3_nplus1;
    p1_prime = p1_pnplus1;
    p2_prime = p2_pnplus1;
    p3_prime = p3_pnplus1;

    

end
    




