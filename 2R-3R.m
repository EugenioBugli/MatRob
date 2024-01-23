clc
format long
syms q1 q2 q3
syms a1 a2 a3 a4
syms l1 l2 l3
syms alpha beta gamma


%% 2R

% direct kine
pos = [l1*cos(q1) + q3*cos(q1+q2); l1*sin(q1) + q3*sin(q1+q2)]
% task variables
p = [pos; q1+q2]
% jacobian
J = jacobian(p, [q1,q2,q3])
simplify(det(J))
Js = simplify(subs(J, [l1, q2], [0.5, pi/2]))
JT = Js.'
tau = JT*[0; 1.5; -4.5]
% position
posin = [2+sqrt(2)/2; sqrt(2)/2]
posfin = [3*sqrt(2)/2; -sqrt(2)/2]
% l1 = 2
% l2 = 1
inverskine(posfin,l1,l2)

%% 3R
clc
% direct kine
pos = [l1*cos(q1) + l2*cos(q1+q2) + l3*cos(q1+q2+q3); l1*sin(q1) + l2*sin(q1+q2) + l3*sin(q1+q2+q3)]
% pos = subs(pos, [l1,l2,l3], [0.5,0.5,0.25])
% pos0 = subs(pos, [q1,q2,q3], [pi/2, -pi/3, 0])
J = jacobian(pos, [q1,q2,q3])
Jsub = subs(J, [q1,q2,q3,l1,l2,l3], [pi/2, -pi/3, 0,0.5,0.5,0.25])
qdot0m = [-pi/6; 0; -pi/2]
v0m = Jsub*qdot0m
double(v0m)
% task variables
% p = [pos; q1+q2]
% % jacobian
% J = simplify(jacobian(p, [q1,q2,q3]))
% simplify(det(J))
% Js = simplify(subs(J, [l1, q2], [0.5, pi/2]))
% JT = Js.'
% tau = JT*[0; 1.5; -4.5]
% position
% % pos = [0.5; sqrt(3)/2]
% % l1 = 0.5
% % l2 = 0.5
% % inverskine(pos,l1,l2)

function inverskine(pos, l1, l2)
    syms q1 q2
    p = [l1*cos(q1) + l2*cos(q1+q2); l1*sin(q1) + l2*sin(q1+q2)];
    c2 = (pos(1)^2 + pos(2)^2 - l1^2 -l2^2)/(2*l1*l2);
    s2_1 = sqrt(1- c2^2);
    s2_2 = -s2_1;
    q2_1 = atan2(s2_1, c2);
    q2_2 = atan2(s2_2, c2);

    det = l1^2 + l2^2 + 2*l1*l2*c2;
    
    s1_1 = (pos(2)*(l1+l2*c2) - pos(1)*l2*s2_1)/det;
    s1_2 = (pos(2)*(l1+l2*c2) - pos(1)*l2*s2_2)/det;

    c1_1 = (pos(1)*(l1+l2*c2) - pos(2)*l2*s2_1)/det;
    c1_2 = (pos(1)*(l1+l2*c2) - pos(2)*l2*s2_2)/det;
   
    q1_1 = atan2(s1_1, c1_1);
    q1_2 = atan2(s1_2, c1_2);

    display('First Solution')

    final1 = [q1_1 q2_1]
    round(rad2deg(final1))
    
    display('Second Solution')
    final2 = [q1_2 q2_2]
    round(rad2deg(final2))
end