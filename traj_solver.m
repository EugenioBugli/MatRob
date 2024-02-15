clc
syms q1 q2 q3 q4 q5 q6
syms a1 a2 a3 a4 a4 a5 a6
syms d0 d1 d2 d3 d4 de l1 l2 l3 l4 N L M N d A B C D K dtcp h p L1 L2
syms alpha beta gamma t

%% quintic example with a0 = 0 and a2 = 0 (3 joints)
pdot0 = [1;-1;0]
T = 2
qg = [0; 0; pi/4]
qs = [-pi/4; pi/4; pi/4]
deltaq = qg-qs

r = [L*cos(q1) + N*cos(q1+q2)*cos(q3);
     L*sin(q1) + N*sin(q1+q2)*cos(q3);
     M + N*sin(q3)]
r = subs(r, [L,M,N], [0.5,0.5,0.5]);
r = jacobian(r, [q1,q2,q3]);
r = subs(r, [q1,q2,q3], qs.');

qdot0 = double(inv(r)*pdot0);

%change this in order to find the specific joint !
fprintf("Solver \n")
for num = 1:3
    if isequal(deltaq(num),0.0)
        display("No changes for joint q"+num)

        p = qs(num)
        subplot(3,1,1)
        fplot(p, LineWidth=1.5)
        title('Position')
        xlabel('s')
        ylabel('m')
        xlim([0 T])
        grid on
        hold on

        pdot = diff(p)
        subplot(3,1,2)
        fplot(0, LineWidth=1.5)

        title('Velocity')
        xlabel('s')
        ylabel('m/s')
        xlim([0 T])
        grid on
        hold on

        pddot = diff(pdot)
        subplot(3,1,3)
        fplot(0, LineWidth=1.5)

        title('Accelleration')
        xlabel('s')
        ylabel('m/s^2')
        xlim([0 T])
        grid on
        hold on
        break

    end
    
    A1 = vpa( (qdot0(num)*T) / deltaq(num) )
    eq1 = a3 + a4 + a5 == 1 - A1;
    eq2 = 3*a3 + 4*a4 + 5*a5 == -A1;
    eq3 = 3*a3 + 6*a4 + 10*a5 == 0;
    [A, b] = equationsToMatrix([eq1;eq2;eq3], [a3,a4,a5]);
    display("Solution for q"+num)
    sol = inv(A)*b

    p = qs(num) + deltaq(num)*(A1*(t/T) + sol(1)*(t/T)^3 + sol(2)*(t/T)^4 + sol(3)*(t/T)^5)
    subplot(3,1,1)
    fplot(p, LineWidth=1.5)

    title('Position')
    xlabel('s')
    ylabel('m')
    xlim([0 T])
    grid on
    hold on

    pdot = diff(p, t)
    subplot(3,1,2)
    fplot(pdot, LineWidth=1.5)

    title('Velocity')
    xlabel('s')
    ylabel('m/s')
    xlim([0 T])
    grid on
    hold on

    pddot = diff(pdot, t)
    subplot(3,1,3)
    fplot(pddot, LineWidth=1.5)

    title('Accelleration')
    xlabel('s')
    ylabel('m/s^2')
    xlim([0 T])
    grid on
    hold on
end

%% cubic example 2R

deltaq = qfin - qin
qi = qin + deltaq*(a*t^3 + b*t^2 + c*t + d)
param = [q1, q2]
m = size(param)
fprintf("Solver \n")
for num = 1:m
    if isequal(deltaq(num),0.0)
        display("No changes for joint q"+num)

        p = qs(num)
        subplot(3,1,1)
        fplot(p, LineWidth=1.5)
        title('Position')
        xlabel('s')
        ylabel('m')
        xlim([0 T])
        grid on
        hold on

        pdot = diff(p)
        subplot(3,1,2)
        fplot(0, LineWidth=1.5)

        title('Velocity')
        xlabel('s')
        ylabel('m/s')
        xlim([0 T])
        grid on
        hold on

        pddot = diff(pdot)
        subplot(3,1,3)
        fplot(0, LineWidth=1.5)

        title('Accelleration')
        xlabel('s')
        ylabel('m/s^2')
        xlim([0 T])
        grid on
        hold on
        break

    end
    
    A1 = vpa( (qdot0(num)*T) / deltaq(num) )
    eq1 = a3 + a4 + a5 == 1 - A1;
    eq2 = 3*a3 + 4*a4 + 5*a5 == -A1;
    eq3 = 3*a3 + 6*a4 + 10*a5 == 0;
    [A, b] = equationsToMatrix([eq1;eq2;eq3], [a3,a4,a5]);
    display("Solution for q"+num)
    sol = inv(A)*b

    p = qs(num) + deltaq(num)*(A1*(t/T) + sol(1)*(t/T)^3 + sol(2)*(t/T)^4 + sol(3)*(t/T)^5)
    subplot(3,1,1)
    fplot(p, LineWidth=1.5)

    title('Position')
    xlabel('s')
    ylabel('m')
    xlim([0 T])
    grid on
    hold on

    pdot = diff(p, t)
    subplot(3,1,2)
    fplot(pdot, LineWidth=1.5)

    title('Velocity')
    xlabel('s')
    ylabel('m/s')
    xlim([0 T])
    grid on
    hold on

    pddot = diff(pdot, t)
    subplot(3,1,3)
    fplot(pddot, LineWidth=1.5)

    title('Accelleration')
    xlabel('s')
    ylabel('m/s^2')
    xlim([0 T])
    grid on
    hold on
end

%% Over-Fly
syms t

A = [3;3]
B = [1;9]
C = [8;9]
v1 = 1
v2 = 2
deltaT = 4

kab = (B-A) / norm(B-A)
kbc = (C-B) / norm(C-B)

Aprime = B - (v1*deltaT*kab)/2

d1 = 0.5*v1*deltaT
d2 = 0.5*v2*deltaT

Cprime = B + d2*kbc

t = linspace(0, deltaT);

poft = Aprime + v1*kab*(t) + 0.5*(t).^2.*(v2*kbc - v1*kab)/deltaT;

% poftdot = diff(poft, t);
% 
% poftddot = diff(poftdot, t);

plot(A(1),A(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
hold on
plot(B(1),B(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
hold on
plot(C(1),C(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
hold on
plot(Aprime(1),Aprime(2),'MarkerEdge',[0.6350 0.0780 0.1840], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
hold on
plot(Cprime(1),Cprime(2),'MarkerEdge',[0.6350 0.0780 0.1840], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
hold on

plot([1 8], [9 9],'k', 'LineWidth', 1.3)
hold on

plot([3 1], [3 9],'k', 'LineWidth', 1.3)
hold on

% plot([Aprime(1) Cprime(1)], [Aprime(2) Cprime(2)],'r', 'LineWidth', 1.3)
% hold on

poft
for i=1:2:199
    elem = poft(i:i+1)
    plot(elem(1), elem(2),'MarkerEdge',[1 0 0], 'LineWidth', 1, 'Marker','.', 'MarkerSize', 18)
end

% subplot(3,1,1)
% fplot(poft, LineWidth=1.5)
% 
% title('Cartesian Path')
% xlabel('s')
% ylabel('(m)')
% xlim([0,deltaT])
% grid on
% 
% subplot(3,1,2)
% fplot(poftdot, LineWidth=1.5)
% 
% title('Cartesian Path')
% xlabel('s')
% ylabel('(m/s)')
% xlim([0,deltaT])
% grid on
% 
% subplot(3,1,3)
% fplot(poftddot, LineWidth=1.5)
% 
% title('Cartesian Path')
% xlabel('s')
% ylabel('(m/s^2)')
% xlim([0,deltaT])
% grid on


%% Rendez-vous
syms t x w
l1 = 0.5
l2 = 0.4
P = [-0.8; 1.1]

plot(P(1),P(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
grid on
hold on
xline(0)
yline(0)
r= l1 + l2 ;
theta=-pi:0.001:pi;
x=r*cos(theta);
y=r*sin(theta);
plot(x,y,"k-", 'LineWidth', 0.6)
val = (-0.8*tand(-20) - 1.1)/tand(-20)
plot([-0.8 val], [1.1 0],'MarkerEdge',[0, 0, 1], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)

k = 0.8*tand(-20) + 1.1
m = tand(-20)
[xo, yo]=linecirc(m, k, 0, 0, 0.9)

intr1 = [xo(1) yo(1)]
intr2 = [xo(2) yo(2)]
plot(intr1(1),intr1(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)
plot(intr2(1),intr2(2),'MarkerEdge',[0.8500 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)

r= 0.6;
theta=-pi:0.001:pi;
x=r*cos(theta)+intr2(1);
y=r*sin(theta)+intr2(2);
plot(x,y,"k-", 'LineWidth', 0.6)

[xo, yo] = linecirc(m,k, intr2(1), intr2(2), 0.6)
intrRV2 = [xo(2) yo(2)]
plot(intrRV2(1),intrRV2(2),'MarkerEdge',[1 0.3250 0.0980], 'LineWidth', 2, 'Marker','.', 'MarkerSize', 25)


%%
syms t T

q1 = pi/4 + (pi/4)*(3*(t/T)^2 - 2*(t/T)^3)
q2 = (-pi/2)*(1 - cos((pi*t)/T))

p1i = subs(q1, t, 0)
p1f = subs(q1, t, T)

p2i = subs(q2, t, 0)
p2f = subs(q2, t, T)

display('Vel')

q1dot = diff(q1,t)
q2dot = diff(q2,t)

p1idot = subs(q1dot, t, 0)
p1fdot = subs(q1dot, t, T)

p2idot = subs(q2dot, t, 0)
p2fdot = subs(q2dot, t, T)

display('Acc')

q1ddot = diff(q1dot,t)
q2ddot = diff(q2dot,t)

p1iddot = subs(q1ddot, t, 0)
p1fddot = subs(q1ddot, t, T)

p2iddot = subs(q2ddot, t, 0)
p2fddot = subs(q2ddot, t, T)

display('Jerk')

q1dddot = diff(q1ddot,t)
q2dddot = diff(q2ddot,t)

