% 
% 
% Author: Eugenio Bugli
% April 2024
% 
% 
clear
clc

syms Q1(t) Q2(t) Q3(t) Q4(t) 
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real
syms l real

% Assume unitary link lenght for this planar 4R

p = l * [cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
    sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];

J = simplify(jacobian(p, [Q1,Q2,Q3]));

time_deriv = [diff(Q1(t),t), diff(Q2(t),t), diff(Q3(t),t), Q1(t), Q2(t), Q3(t)];

% take the time derivative of the jacobian

oper = diff(J, t);

% substitute symbolic variable to reduce expressions

dJ = expand(subs(oper, time_deriv , [dq1,dq2,dq3, q1,q2,q3]));

% apply some collect, expand and simplify to have a simpler expression
% you may have to change the parameters to collect

dJdq = collect(simplify(dJ * [dq1,dq2,dq3].'), [dq1,dq2,dq3])