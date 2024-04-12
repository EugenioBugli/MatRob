clc
syms q1 q2 q3 l d real
syms dq1 dq2 dq3 real

p = [l*cos(q1)+ l*cos(q1+q2) + l*cos(q1+q2+q3); l*sin(q1) + l*sin(q1+q2) + l*sin(q1+q2+q3)];

Jrcm = jacobian(p, [q1,q2,q3]) % 2x3

Jact = Jrcm(:,1) % 2x1
Junact = [-Jrcm(:,2) -Jrcm(:,3)] %2x2


% Jact * dq1 = Junact * [dq2; dq3]


% directly without substitution
simplify(det(Junact))  

%  with substitution
calc = simplify(subs(Junact, [q2,q3], [-q1, pi - q1]));
simplify(det(calc));

simplify(inv(Junact))

theta_unact_dot = simplify(inv(Junact)*Jact)

theta3_dot = theta_unact_dot(2)

simplify(subs(theta3_dot, [q2,q3], [-q1, q1 - pi ]))