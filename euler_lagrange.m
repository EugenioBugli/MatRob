% 
% 
% Author: Eugenio Bugli
% April 2024
% 
% 
clear
clc
syms m1 m2 m3 m4 real
syms q1 q2 q3 q4 real
syms a1 a2 a3 a4 a5 a6 real
syms q1dot q2dot q3dot q4dot real

syms dc1 dc2 real
syms l real
syms Ic1 Ic2 real
syms g0 real

T1 = 0.5*(m1*q1^2 + Ic1)*q1dot^2;
vc2 = [q1dot - l*q2dot*sin(q2); q2dot*l*cos(q2)];
T2 = 0.5*m2*sq_norm(vc2) + 0.5*(Ic2 + m2*dc2^2)*q2dot^2;

T = simplify(expand(T1+T2));

qdot = [q1dot; q2dot];

M = getInertiaMatrix(T, qdot)

qdot = [q1dot; q2dot];
q = [q1, q2];
c = getCorio_Centrif(M, qdot, q)

U = -m1*g0*q1 -m2*g0*(q1+dc2*cos(q2));

g = jacobian(U,q)

function c = getCorio_Centrif(M, qd, q)
    disp("Coriolis terms: i!=j");
    disp("Centrifugal terms: i=j");
    c = simplify(getChristoffel(M, qd, q));
end

function C = getChristoffel(M, qd, q)
    s = size(M,2);
    C = sym(zeros(s,1));
    for i=1:s
        Mkdiff = jacobian(M(:,i), q);
        res = 0.5 * qd.' * ( Mkdiff + Mkdiff.' - diff(M, q(i)) ) * qd;
        C(i) = simplify(res);
    end
end

function M = getInertiaMatrix(T, q)
    s = size(q,1);
    M = sym(zeros(s));
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

function val = sq_norm(vec)
    n = norm(vec)^2;
    val = simplify(n);
end