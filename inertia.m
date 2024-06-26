clear
clc
syms m1 m2 m3 m4 real
syms q1 q2 q3 q4 real
syms a1 a2 a3 a4 a5 a6 real
syms q1dot q2dot q3dot q4dot real


vc1 = [q1dot; 0];
vc2 = [q1dot; q2dot];
vc3 = [q1dot + q3dot; q2dot];
vc4 = [q1dot + q3dot; q2dot + q4dot];


T1 = 0.5*m1*sq_norm(vc1)
T2 = 0.5*m2*sq_norm(vc2)
T3 = 0.5*m3*sq_norm(vc3)
T4 = 0.5*m4*sq_norm(vc4)

T = simplify(expand(T1+T2+T3+T4))

% qdot = [q1dot; q2dot; q3dot; q4dot];

% M = getInertiaMatrix(T, qdot)

qdot = [q1dot; q2dot; q3dot];
q = [q1, q2, q3];

M = [a1 + 2*a2*q2 + a3*q2^2 + 2*a4*q2*sin(q3) + a5*(sin(q3))^2, 0, 0;
    0, a3, a4*cos(q3);
    0, a4*cos(q3), a6]

c = getCorio_Centrif(M, qdot, q)

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