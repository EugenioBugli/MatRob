clc
format long
syms q1 q2 q3 q4
syms a1 a2 a3 a4
syms l1 l2 l3
syms alpha beta gamma

pos = [l1 - l2*cos(q1); -l2*sin(q1); 2*q1 + q2]
posmin = [l1 - l2*cos(q1); -l2*sin(q1); pi + q1]

J = jacobian(pos, [q1,q2,q3,q4])

rank(J)

Jmin = jacobian(posmin, [q1,q2,q3,q4])

rank(Jmin)

%% 2R

% direct kine
pos = [l1*cos(q1) + l2*cos(q1+q2); l1*sin(q1) + l2*sin(q1+q2)] 
% jacobian
% po = subs(pos, [l1,l2], [1,1])
% J = jacobian(po, [q1,q2])


posfin = [5; 1]
l1 = 2
l2 = 1

inverskine_2(posfin,l1,l2)

%% 3R
clc
% direct kine
r = [l1*cos(q1) + l2*cos(q1+q2) + l3*cos(q1+q2+q3); l1*sin(q1) + l2*sin(q1+q2) + l3*sin(q1+q2+q3)]%; q1+q2+q3]
r = l1*[cos(q1) + cos(q1+q2) + cos(q1+q2+q3); sin(q1) + sin(q1+q2) + sin(q1+q2+q3)]
% pos = subs(pos, [l1,l2,l3], [0.5,0.5,0.25])
% pos0 = subs(pos, [q1,q2,q3], [pi/2, -pi/3, 0])
J = jacobian(r, [q1,q2,q3])
Jsub = subs(J, [q1,q2,q3,l1], [0,0,pi/2,0.5])
pinv(Jsub)
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
% % inverskine_3(pos,l1,l2,l3)

function inverskine_2(pos, l1, l2)
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

function inverskine_3(pos, l1, l2, l3)
    syms q1 q2 q3
    r = [l1*cos(q1) + l2*cos(q1+q2) + l3*cos(q1+q2+q3); l1*sin(q1) + l2*sin(q1+q2) + l3*sin(q1+q2+q3); q1+q2+q3];

    c2 = (pos(1)^2 + pos(2)^2 + l3^2 - 2*l3*( pos(1)*cos(pos(3)) + pos(2)*sin(pos(3)) ) - l1^2 - l2^2 )/(2*l1*l2);
    s2_1 = sqrt(1- c2^2);
    s2_2 = -s2_1;
    q2_1 = atan2(s2_1, c2);
    q2_2 = atan2(s2_2, c2);

    s1_1 = (l1 + l2*c2)*(pos(2) - l3*sin(pos(3))) + l2*s2_1*(pos(1) - l3*cos(pos(3)));
    s1_2 = (l1 + l2*c2)*(pos(2) - l3*sin(pos(3))) + l2*s2_2*(pos(1) - l3*cos(pos(3)));

    c1_1 = (l1 + l2*c2)*(pos(1) - l3*cos(pos(3))) + l2*s2_1*(pos(2) - l3*sin(pos(3)));
    c1_2 = (l1 + l2*c2)*(pos(1) - l3*cos(pos(3))) + l2*s2_2*(pos(2) - l3*sin(pos(3)));
   
    q1_1 = atan2(s1_1, c1_1);
    q1_2 = atan2(s1_2, c1_2);

    q3_1 = pos(3) - q1_1 - q2_1;
    q3_2 = pos(3) - q1_2 - q2_2;

    display('First Solution')

    final1 = [q1_1 q2_1 q3_1]
    round(rad2deg(final1))
    
    display('Second Solution')
    final2 = [q1_2 q2_2 q3_2]
    round(rad2deg(final2))
end