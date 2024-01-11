% Ri = [sqrt(2)/2,0,sqrt(2)/2;
%       0,-1,0;
%       sqrt(2)/2,0,-sqrt(2)/2]
% Rf = x_m(-pi/2)*y_m(pi/3)*z_m(pi/3)
% 
% rot = Ri.'*Rf
% 
% check(rot)
format long
syms q1 q2 q3 q4 q5 q6
syms a1 a2 a3 a4
syms d0 d1 d2 d3 d4 de l1 l2 l3 l4 N L M N d A B C D K dtcp h p L1 L2 t

Jac = [-sin(q1)*(l2*cos(q2)+l3*cos(q3)) -l2*cos(q1)*sin(q2) -l3*cos(q1)*sin(q3);
      cos(q1)*(l2*cos(q2)+l3*cos(q3)) -l2*sin(q1)*sin(q2) -l3*sin(q1)*sin(q3);
      0 l2*cos(q2) l3*cos(q3)]

R01 = z_m(q1)
J1 = simplify(R01.' * Jac)
simplify(simplify(det(J1)))

rank2 = simplify(subs(Jac, [q2], [q3]))
rank(rank2)

rank1 = simplify(subs(rank2, [q3], [pi/2]))
rank(rank1)

null(rank1)
colspace(rank1)
null(rank1.')
colspace(rank1.')

% % % T01 = direct(sym([q1]),sym([pi/2]),[A],[B],1)
% % % R01 = T01(1:3,1:3)
% % % J1 = simplify(R01.' * Jacobian)
% % % subs(J1, [B,D], [0,0])

% rangespace = simplify(colspace(sing))
% nullspace = simplify(null(sing))


% % sing = simplify(subs(Jacobian, [B,D], [0,0]))
% % Jac1 = simplify(R01.'*sing)
% % m = simplify(sing*sing.')
% % simplify(det(m))
% % 
% % nullspace = simplify(null(sing))
% % rangespace = simplify(colspace(sing))
% % nulltransposed = simplify(null(sing.'))

% R04 = direct(sym([0,q2,0,q4]),sym([0,-pi/2,0,0]),[q1,0,q3,0],[0,N,0,0],4)
% p04 = R04(1:3,4)


% % T02 = direct(sym([q1,q2]),sym([0,-pi/2]),[0,0],[L1,L2],2)
% % simplify(T02)
% % 
% % T02n = direct(sym([pi/2,-pi/2]),sym([0,-pi/2]),[0,0],[0.5,0.6],2)

% T02 = direct([0,q2,q3],sym([pi/2,0,0]),[q1,0,0],[0,L,L],3);
% simplify(T02);

% pos = [q1 + L*cos(q2) + L*cos(q2+q3); L*sin(q2) + L*sin(q2+q3); q2+q3]
% Jacobian = [diff(pos,q1),diff(pos,q2), diff(pos,q3)]
% jacobian(pos, [q1,q2,q3])
% pos = [q2*cos(q1) + q4*cos(q1+q3); q2*sin(q1) + q4*sin(q1+q3); q1+q3]
% Jacobian = jacobian(pos, [q1,q2,q3,q4])
% det(Jacobian)
% simplify(det(simplify(Jacobian*Jacobian.')))
% sing = simplify(subs(Jacobian, [q2,q3], [0,0]))
% sing = Jacobian
% rank = rank(sing)
% nullspace = simplify(null(sing))
% rangespace = simplify(colspace(sing))
% nulltransposed = simplify(null(sing.'))
% 
% r = [1;0;0]
% inv = simplify(pinv(sing))
% 
% q = inv * r
% 
% ver = simplify(sing*q)

% v = [L*cos(q1); K*q1; L*sin(q1)]
% simplify(norm(v))

function mat = direct(theta,alpha,d,a,num)

    mat = dh(alpha(1),theta(1),d(1),a(1))

    for i = 2:num
        dh(alpha(i),theta(i),d(i),a(i))
        mat = mat * dh(alpha(i),theta(i),d(i),a(i));
        simplify(mat);
    end
end

function mdh = dh(alpha,t,d,a)
    mdh = [cos(t),  -cos(alpha)*sin(t),     sin(alpha)*sin(t),     a*cos(t);
           sin(t),   cos(alpha)*cos(t),    -sin(alpha)*cos(t),     a*sin(t);
           0,               sin(alpha),            cos(alpha),            d;
           0                     0,                 0,            1];
end

function bool = check(m)
    deter = int64(det(m)) == 1;
    nor1 = int64(norm(m(:,1))) == 1;
    nor2 = int64(norm(m(:,2))) == 1;
    nor3 = int64(norm(m(:,3))) == 1;
    orthog = isequal(inv(m),m.');
    bool = deter & nor1 & nor2 & nor3 & orthog;
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

function s = skew(r)
s = [0,-r(3),r(2);
     r(3),0,-r(1);
     -r(2),r(1),0];
end