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









