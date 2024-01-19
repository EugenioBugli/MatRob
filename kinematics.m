% Ri = [sqrt(2)/2,0,sqrt(2)/2;
%       0,-1,0;
%       sqrt(2)/2,0,-sqrt(2)/2]
% Rf = x_m(-pi/2)*y_m(pi/3)*z_m(pi/3)
% 
% rot = Ri.'*Rf
% 
% check(rot)
clc
format long
syms q1 q2 q3 q4 q5 q6
syms a1 a2 a3 a4
syms d0 d1 d2 d3 d4 de l1 l2 l3 l4 N L M N d A B C D K dtcp h p L1 L2
syms t q1(t) q2(t) q3(t) alpha beta gamma

R = y_m(q3(t))*z_m(q2(t))*z_m(q1(t))
Rdot = simplify(diff(R, t))
Rdot1 = simplify(subs(Rdot, [diff(q1(t),t),diff(q2(t),t),diff(q3(t),t)], [alpha,beta,gamma]))
Rdot2 = simplify(subs(Rdot1, [q1(t),q2(t),q3(t)], [d1, d2, d3]))
R2 = simplify(subs(R, [q1(t),q2(t),q3(t)], [d1, d2, d3]))
skews = simplify(Rdot2* R2.')
w = skewTovec(skews)

x = [alpha, beta, gamma]
A = equationsToMatrix(w, x)

ObtainfromOrientation(R, 'YZX', [t, q1(t), q2(t), q3(t), alpha, beta, gamma, d1, d2, d3], 'RPY')

% Jac = [-sin(q1)*(l2*cos(q2)+l3*cos(q3)) -l2*cos(q1)*sin(q2) -l3*cos(q1)*sin(q3);
%       cos(q1)*(l2*cos(q2)+l3*cos(q3)) -l2*sin(q1)*sin(q2) -l3*sin(q1)*sin(q3);
%       0 l2*cos(q2) l3*cos(q3)]
% 
% R01 = z_m(q1)
% J1 = simplify(R01.' * Jac)
% simplify(simplify(det(J1)))
% 
% rank2 = simplify(subs(Jac, [q2], [q3]))
% rank(rank2)
% 
% rank1 = simplify(subs(rank2, [q3], [pi/2]))
% rank(rank1)
% 
% null(rank1)
% colspace(rank1)
% null(rank1.')
% colspace(rank1.')

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
function ObtainfromOrientation(R, seq, vars, typ)
    display('If RPY double check the correct order!')
    t = vars(1);
    q2(t) = vars(3);
    if typ == 'RPY'
        fprintf('RPY case so: \n q1(t) is gamma \n q2(t) is beta \n q3(t) is alpha \n')
        q1(t) = vars(4);
        q3(t) = vars(2);
        alpha = vars(7);
        gamma = vars(5);
        d1 = vars(10);
        d3 = vars(8);
    else
        q1(t) = vars(2);
        q3(t) = vars(4);
        alpha = vars(5);
        gamma = vars(7);
        d1 = vars(8);
        d3 = vars(10);
    end
    
    beta = vars(6);
    d2 = vars(9);

    if seq == 'XYZ'
        m = x_m(q1(t))*y_m(q2(t))*z_m(q3(t))
    elseif seq == 'XZY'
        m = x_m(q1(t))*z_m(q2(t))*y_m(q3(t))
    elseif seq == 'XYX'
        m = x_m(q1(t))*y_m(q2(t))*x_m(q3(t))
    elseif seq == 'XZX'
        m = x_m(q1(t))*z_m(q2(t))*x_m(q3(t))

    elseif seq == 'YXZ'
        m = y_m(q1(t))*x_m(q2(t))*z_m(q3(t))
    elseif seq == 'YZX'
        m = y_m(q1(t))*z_m(q2(t))*x_m(q3(t))
    elseif seq == 'YXY'
        m = y_m(q1(t))*x_m(q2(t))*y_m(q3(t))
    elseif seq == 'YZY'
        m = y_m(q1(t))*z_m(q2(t))*y_m(q3(t))

    elseif seq == 'ZXY'
        m = z_m(q1(t))*x_m(q2(t))*y_m(q3(t))
    elseif seq == 'ZYX'
        m = z_m(q1(t))*y_m(q2(t))*x_m(q3(t))
    elseif seq == 'ZXZ'
        m = z_m(q1(t))*x_m(q2(t))*z_m(q3(t))
    elseif seq == 'ZYZ'
        m = z_m(q1(t))*y_m(q2(t))*z_m(q3(t))
    end
    
    mdot = simplify(diff(m,t))
    display('derivatives of q1(t), q2(t), q3(t) are now called alpha, beta, gamma')
    mdot1 = simplify(subs(mdot, [diff(q1(t),t),diff(q2(t),t),diff(q3(t),t)], [alpha,beta,gamma]));
    display('q1(t), q2(t), q3(t) are the now d1, d2, d3')
    mdot2 = simplify(subs(mdot1, [q1(t),q2(t),q3(t)], [d1, d2, d3]));
    
    if typ == 'RPY'
        fprintf('RPY case so: \n q1dot(t) is gamma \n q2dot(t) is beta \n q3dot(t) is alpha \n\n q1(t) is d3 \n q2(t) is d2 \n q3(t) is d1 \n')
    end
    
    m2 = simplify(subs(m, [q1(t),q2(t),q3(t)], [d1, d2, d3]))
    display('Compute S(w) = Rdot * R^T \n')
    S = simplify(mdot2*m2.')
    display('decompose it in a vector \n')
    w = skewTovec(S)
    display('Given Ax = b with x = [alpha, beta, gamma] and b = w, here you have A \n')
    if typ == 'RPY'
        x = [gamma, beta, alpha]
    else 
        x = [alpha, beta, gamma]
    end

    A = equationsToMatrix(w, x)
    
    display('find singularities \n')
    d = simplify(det(A))
    solution = solve(d, 0)
    
    display('A with singularities')
    Asing = subs(A, d2, pi/2)

    display('angular velocity not realizable (not belong to R(tphi))')
    null(Asing.')
    display('time derivatives that generate zero angular velocity (belongs to N(tphi))')
    null(Asing)

end

function mat = GeometricJacobian(pos, param)
end

function mat = AnalyticJacobian(pos, param)
    display('Normal Jacobian')
    mat = jacobian(pos, param)
    simplify(det(mat))
    R01 = z_m(q1);
    display('Jacobian expressed in frame 1')
    mat1 = simplify(R01.' * mat)
    simplify(det(mat1))
    solve(simplify(det(mat1))==0)
end

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

function w = skewTovec(S)
w = [S(3,2); S(1,3);S(2,1)];
end