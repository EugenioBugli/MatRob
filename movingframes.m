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

% use diag matrices for each joint then multiply

Ic1 = diag([Ic1xx, Ic1yy, Ic1zz]);
Ic2 = diag([Ic2xx, Ic2yy, Ic2zz]);
Ic3 = diag([Ic3xx, Ic3yy, Ic3zz]);
Ic4 = diag([Ic4xx, Ic4yy, Ic4zz]);

alpha = [0,0];
a = [l1, l2];
d = [0, 0];
theta = [q1, q2];

DH = [sym(alpha); a; d; sym(theta)];

sig = [0,0]; % fill these values with 0 if your joint is Prismatic, otherwise fill with 0
qdot = [dq1,dq2];
mass = [m1,m2];
iner = [Ic1 Ic2];

COM_Pos = [[-l1 + d1;0;0],[-l2 + d2;0;0]];

VEL = comp_velocities(DH, sig, qdot);

AngVel = VEL(1:3,1:2)
LinVel = VEL(1:3,3:4)

COM_Vel = comp_com_vel(AngVel, LinVel, COM_Pos)

T = comp_kine(AngVel, COM_Vel, mass, iner)

M = simplify(getInertiaMatrix(T, qdot))

function CVel = comp_com_vel(ang, lin, pos)
    num = size(ang,2);
    CVel = [];
    for i=1:num
        vci = lin(:,i) + cross(ang(:,i),pos(:,i));
        CVel  = [CVel vci];
    end
    CVel = simplify(CVel);
end

function T = comp_kine(ang, lin, mass, iner)
    num = size(ang,2);
    T = 0;
    for i=1:num
        lin_part = 0.5*mass(i)*sq_norm(lin(:,i));
        ang_part = 0.5*ang(:,i).'*iner(:,(1 + (3*(i-1))):(3 + (3*(i-1))))*ang(:,i);
        Ti = simplify(lin_part + ang_part);
        T = simplify(T + lin_part + ang_part);
    end
end

function OUT = comp_velocities(DH, sig, qdot)
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
    v_00 = 0;

    w_11 = mat(1:3,1:3).' * (w_00 +  (1-sig(1)) * [0;0;qdot(1)]);
    ang = [w_11];

    v_11 = mat(1:3,1:3).' * (v_00 +  sig(1)*[0;0;qdot(1)]) + cross(w_11,  mat(1:3,1:3).' * mat(1:3,4));
    lin = [v_11];
    
    for i = 2:num
        % i-1^A_i
        rot = dh(alpha(i),theta(i),d(i),a(i));

        % display("Angular Velocity W_"+i+"_"+i);
        w = rot(1:3,1:3).' * ( ang(:,i-1) +  (1-sig(i)) * [0;0;qdot(i)]);

        % display("Linear Velocity V_"+i+"_"+i);
        v = rot(1:3,1:3).' * ( lin(:,i-1) +  sig(i) * [0;0;qdot(i)]) + cross(w,  rot(1:3,1:3).' * rot(1:3,4));

        ang = [ang w];
        lin = [lin simplify(v)];

        % 0^A_i
        mat = simplify(mat * dh(alpha(i),theta(i),d(i),a(i)));
    end
    OUT = simplify([ang, lin]);
end

function mdh = dh(alpha,t,d,a)
    mdh = [cos(t),  -cos(alpha)*sin(t),     sin(alpha)*sin(t),     a*cos(t);
           sin(t),   cos(alpha)*cos(t),    -sin(alpha)*cos(t),     a*sin(t);
           0,               sin(alpha),            cos(alpha),            d;
           0                     0,                 0,            1];
end

function val = sq_norm(vec)
    n = norm(vec)^2;
    val = simplify(n);
end

function M = getInertiaMatrix(T, q)
    s = size(q,2)
    for i=1:s
        for j=i:s
            if i == j
                M(i,j) = diff(T,q(i),2);
            else
                temp = diff(T,q(i));
                M(i,j) = diff(temp, q(j));
                M(j,i) = M(i,j); % symmetry
            end
        end
    end
end

% function mat = direct(DH)
%     alpha = DH(1,:);
%     a = DH(2,:);
%     d = DH(3,:);
%     theta = DH(4,:);
% 
%     num = size(theta,2);
%     numbers = [1:num]
%     table_Dh = table(numbers.', alpha.', a.', d.', theta.','VariableNames', ["i", "alpha", "a", "d", "theta"])
% 
%     display("Matrix T01")
% 
%     mat = dh(alpha(1),theta(1),d(1),a(1))
% 
%     for i = 2:num
%         display("Matrix T"+int2str(i-1)+":"+i)
%         dh(alpha(i),theta(i),d(i),a(i))
%         display("Matrix T0"+i)
%         mat = simplify(mat * dh(alpha(i),theta(i),d(i),a(i)))
%     end
% end