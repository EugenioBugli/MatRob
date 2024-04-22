% 
% 
% Author: Eugenio Bugli
% April 2024
% 
% 

clear
clc

syms q1 q2 q3 q4 real
syms a1 a2 a3 a4 real
syms d1 d2 d3 d4 real
syms l1 l2 l3 l4 real
syms dq1 dq2 dq3 dq4 real
syms m1 m2 m3 m4 real
syms Ic1xx Ic2xx Ic3xx Ic4xx real
syms Ic1yy Ic2yy Ic3yy Ic4yy real
syms Ic1zz Ic2zz Ic3zz Ic4zz real

% syms Ic1 Ic2 Ic3 Ic4 real % comment here if you need to use matrix form for Ici


alpha = [0,pi/2,0];
a = [l1, 0, l3];
d = [0, l2, 0];
theta = [q1, q2, q3];

DH = [sym(alpha); a; d; sym(theta)];

sig = [0,0, 0]; % fill these values with 1 if your joint is Prismatic, otherwise fill with 0
qdot = [dq1,dq2, dq3];
mass = [m1,m2, m3];
iner = [Ic1 Ic2 Ic3];

N = 3

% CoM position w.r.t origin of frame i:
COM_Pos = [[-l1 + d1; 0; 0],[0; -l2 + d2; 0], [-l3 + d3; 0; 0]];

% FORWARD RECURSION
OUT_FR = simplify(forward_recursion(DH, sig, qdot, COM_Pos));

AngVel = OUT_FR(1:3,1:N)
AngAcc = OUT_FR(1:3,(N+1):2*N)
LinAcc = OUT_FR(1:3,(2*N+1):3*N)
CoMAcc = OUT_FR(1:3,(3*N+1):4*N)

% BACKWARD RECURSION
OUT_BR = backward_recursion(DH, sig, M, AngVel, AngAcc, CoMAcc, I)

Forces = OUT_FR(1:3,1:N)
Moments = OUT_BR(1:3, (N+1):2*N)

% FORCE PROJECTION
OUT_FP = force_projection(DH, Forces, Moments, Dissp)

function OUT = force_projection(DH, sig, Forces, Moments, Dissp, qdot)
    alpha = DH(1,:);
    a = DH(2,:);
    d = DH(3,:);
    theta = DH(4,:);
    
    num = size(theta,2);
    numbers = [1:num];
    table_Dh = table(numbers.', alpha.', a.', d.', theta.','VariableNames', ["i", "alpha", "a", "d", "theta"])

    % display("Matrix T01");

    mat = dh(alpha(1),theta(1),d(1),a(1));

    gen_forces = sym(zeros(num));
    
    for i = 1:num
        % i^A_i+1
        rot = dh(alpha(i),theta(i),d(i),a(i));

        % display("Gen_Force U_"+i+"_"+i);
        
        % you may need to change this according to the dissipative terms
        % that you have. In this case only Viscous Friction is used.
        if sig(i)==1
            % Prismatic Joint
            u = Forces(i).' * [0;0;1] + Dissp(i)*qdot(i)
        else
            % Revolute Joint
            u = Moments(i).' * [0;0;1] + Dissp(i)*qdot(i)
        end

        % i^z_i-1 = (i-1^R_i)^T * i-1^z_i-1 = [0;0;1]
        
        gen_forces(i) = u;

        % 0^A_i
        mat = simplify(mat * dh(alpha(i),theta(i),d(i),a(i)));
    end
    OUT = simplify(gen_forces);
end

function OUT = backward_recursion(DH, sig, M, AngVel, AngAcc, CoMAcc, I)
    alpha = DH(1,:);
    a = DH(2,:);
    d = DH(3,:);
    theta = DH(4,:);
    
    num = size(theta,2);
    numbers = [1:num];
    table_Dh = table(numbers.', alpha.', a.', d.', theta.','VariableNames', ["i", "alpha", "a", "d", "theta"])

    % display("Matrix T01");

    mat = dh(alpha(1),theta(1),d(1),a(1));

    forc = sym(zeros(num));
    mome = sym(zeros(num));
    
    for i = num:1
        % i^A_i+1
        rot = dh(alpha(i),theta(i),d(i),a(i));

        % display("Force F_"+i+"_"+i);
        f = rot(1:3,1:3).'*forc(num-i+1) + M(i)*(CoMAcc(i))

        % display("Moment Tau_"+i+"_"+i);
        tau = rot(1:3,1:3).'*mome(num-i+1) + cross( (rot(1:3,1:3).'*forc(num-i+1)), pos(i)) - cross( f, (rot(1:3,4)+pos(i)) ) + I(i)*AngAcc(i) + cross( AngVel(i), (I(i)*AngVel(i)))

        forc(i) = f;
        mome(i) = tau;

        % 0^A_i
        mat = simplify(mat * dh(alpha(i),theta(i),d(i),a(i)));
    end
    OUT = simplify([forc, mome]);
end

function OUT = forward_recursion(DH, sig, qdot, pos)
    alpha = DH(1,:);
    a = DH(2,:);
    d = DH(3,:);
    theta = DH(4,:);
    
    num = size(theta,2);
    numbers = [1:num];
    table_Dh = table(numbers.', alpha.', a.', d.', theta.','VariableNames', ["i", "alpha", "a", "d", "theta"])

    % display("Matrix T01");

    mat = dh(alpha(1),theta(1),d(1),a(1));

    w_00 = 0;
    dw_00 = 0;
    a_00 = 0 - g0;

    w_11 = mat(1:3,1:3).' * (w_00 + [0;0;qdot(1)]);
    ang = [w_11];

    dw_11 = mat(1:3,1:3).' * (dw_00 + [0;0;qddot(1)] + qdot(1)*cross(w_00,[0;0;1]));
    ang_der = [dw_11];

    a_11 = mat(1:3,1:3).' * a_00 + cross(dw_11,mat(1:3,4)) + cross(w_11, cross(w_11, mat(1:3,4)));
    acc = [a_11];

    ac_11 = a_11 + cross( dw_11,pos(1) ) + cross(w_11, cross(w_11, pos(1)));
    acc_com = [ac_11];
    
    for i = 2:num
        % i-1^A_i
        rot = dh(alpha(i),theta(i),d(i),a(i));

        % display("Angular Velocity W_"+i+"_"+i);
        w = rot(1:3,1:3).' * ( ang(:,i-1) + [0;0;qdot(i)]);

        % display("Angular Acceleration dW_"+i+"_"+i);
        dw = rot(1:3,1:3).' * ( ang_der(:,i-1) + [0;0;qddot(i)] + qdot(i)*cross(ang(:,i-1),[0;0;1]));

        % display("Linear Acceleration A_"+i+"_"+i);
        a = rot(1:3,1:3).' * acc(:,i-1) + cross(dw,rot(1:3,4)) + cross(w, cross(w, rot(1:3,4)));

        % display("CoM's Acceleration AC_"+i+"_"+i);
        acom = a + cross( dw,pos(i) ) + cross(w, cross(w, pos(i)));

        ang = [ang w];
        ang_der = [ang_der simplify(dv)];
        acc = [acc simplify(a)];
        acc_com = [acc_com simplify(acom)];

        % 0^A_i
        mat = simplify(mat * dh(alpha(i),theta(i),d(i),a(i)));
    end
    OUT = simplify([ang, ang_der, acc, acc_com]);
end

function mdh = dh(alpha,t,d,a)
    mdh = [cos(t),  -cos(alpha)*sin(t),     sin(alpha)*sin(t),     a*cos(t);
           sin(t),   cos(alpha)*cos(t),    -sin(alpha)*cos(t),     a*sin(t);
           0,               sin(alpha),            cos(alpha),            d;
           0                     0,                 0,            1];
end