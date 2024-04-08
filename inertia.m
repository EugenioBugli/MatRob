clear all
clc
syms m1 m2 m3 m4 real
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

qdot = [q1dot; q2dot; q3dot; q4dot];

M = getInertiaMatrix(T, qdot)

function M = getInertiaMatrix(T, q)
    s = size(q,1);
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