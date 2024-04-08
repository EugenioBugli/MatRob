clc
syms q1 q2 q3 q4 q5 q6 q7
syms L L0 L1 L2 L3 a b
syms d1 d2 d3 d4 d5 d6 d7
syms a1 a2 a3 a4

T03 = simplify(direct(sym([0,pi/2,q3]), sym([pi/2,pi/2,0]),[q1,q2,L], [0,0,0], 3))

r = T03(1:3,4)

rsub = subs(r, [q1,q2,q3], [0, pi/6, -pi/2])
vpa(rsub)

JL = simplify(jacobian(r, [q1,q2,q3]))
det(JL)
Js = subs(JL, q3, 0)
null(Js)
null(Js.')
% JA = [0 sin(q1) sin(q1) sin(q2+q3)*cos(q1); 0 -cos(q1) -cos(q1) sin(q2+q3)*sin(q1); 1 0 0 -cos(q2+q3)]
% 
% Jg = [JL; JA]
% 
% R01 = [cos(q1), 0,  sin(q1); sin(q1), 0, -cos(q1); 0, 1,        0];
% 
% J1L = simplify(R01.' * JL)
% 
% rank(J1L)
% 
% simplify(det(J1L(1:2,2:3)))
% 
% Jsubs = subs(J1L, [q2,q3], [0,0])
% rank(Jsubs)
% 
% J1A = simplify(R01.' * JA);
% 
% J1g = [J1L; J1A];

function mat = direct(theta,alpha,d,a,num)

    numbers = [1:num];
    table_Dh = table(numbers.', alpha.', a.', d.', theta.','VariableNames', ["i", "alpha", "a", "d", "theta"])
    display("Matrix T01")
    mat = dh(alpha(1),theta(1),d(1),a(1))
    mat(1:3,3)
    for i = 2:num
        display("Matrix T"+int2str(i-1)+":"+i)
        dh(alpha(i),theta(i),d(i),a(i))
        display("Matrix T0"+i)
        mat = mat * dh(alpha(i),theta(i),d(i),a(i))
        simplify(mat);
        mat(1:3,3)
    end
end

function mdh = dh(alpha,t,d,a)
    mdh = [cos(t),  -cos(alpha)*sin(t),     sin(alpha)*sin(t),     a*cos(t);
           sin(t),   cos(alpha)*cos(t),    -sin(alpha)*cos(t),     a*sin(t);
           0,               sin(alpha),            cos(alpha),            d;
           0                     0,                 0,            1];
end

function mat = z_m(a)
mat = [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1];
end
function mat = y_m(a) 
mat = [cos(a), 0, sin(a); 0, 1, 0; -sin(a), 0, cos(a)];
end

function mat = x_m(a) 
mat = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
end