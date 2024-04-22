% 
% 
% Author: Eugenio Bugli
% April 2024
% 
% 
clear
clc

syms q1 q2 real
syms dc2 real
syms dq1 dq2 real
syms Ic2 real
syms m1 m2 real


% example of factorizations for c(q,dq):

M = [m1 + m2, - m2*dc2*sin(q2); - m2*dc2*sin(q2), Ic2 + m2*dc2^2];
dM = [0, - dq2*m2*dc2*cos(q2); - dq2*m2*dc2*cos(q2), 0]

C1 = [0 0; 0 -m2*dc2*cos(q2)]
C2 = zeros(2)

c1 = [dq1,dq2]*[C1]*[dq1;dq2];
c2 = [dq1,dq2]*[C2]*[dq1;dq2];

c = [c1;c2]

%% Use Christoffel Symbols

% compute the i-th rows of S(q,dq) by doing the product :
% dq^T * C_i(q)  with dq^T you mean a row vector!

S = [[dq1,dq2]*C1;[dq1,dq2]*C2]

% dM - 2*S is a skew-symmetric :

skew_chr = dM - 2*S

% skew symmetry property check :

-skew_chr == skew_chr.'

% dq^T *(dM - 2S_prime) * dq = 0 check:

simplify([dq1,dq2]*skew_chr*[dq1;dq2])
%% By construction
syms Y real
% I am adding this quantities because by doing S_prime * [dq1;dq2] we get
% again c(q, dq)

S_prime = S + [0, 0; dq2 -dq1]

% you can do also multiply something (Y) to the matrix added to S.
% S_prime = S + Y*[0, 0; dq2 -dq1] and will work 
% if S was a correct factorization of c(q,dq) then also S_prime will be correct

S_prime*[dq1;dq2]

skew_prime = dM - 2*S_prime

% NOT a skew-symmetric matrix

-skew_prime == skew_prime.'

% However, it is a valid factorization because satisfies this:
% dq^T *(dM - 2S_prime) * dq = 0

simplify([dq1,dq2]*skew_prime*[dq1;dq2])


%% Example from Midterm 23:
clear
clc

syms q1 q2 q3 real
syms a1 a2 a3 a4 a5 a6 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real

M = [a1 + 2*a2*q2 + a3*q2^2 + 2*a4*q2*sin(q3) + a5*sin(q3)^2 0 0;
     0 a3 a4*cos(q3); 0 a4*cos(q3) a6]

dM = [2*a2*dq2 + 2*a3*q2*dq2 + 2*a4*dq2*sin(q3) + 2*a4*q2*cos(q3)*dq3 + 2*a5*sin(q3)*cos(q3)*dq3 0 0; 0 0 -a4*dq3*sin(q3); 0 -a4*dq3*sin(q3) 0]

c = simplify(getCorio_Centrif(M, [dq1;dq2;dq3], [q1,q2,q3]));

s = size(M,2);
C = sym(zeros(3*s,s));
q = [q1,q2,q3];
for i=1:s
    Mkdiff = jacobian(M(:,i), q);
    res = 0.5 * ( Mkdiff + Mkdiff.' - diff(M, q(i)) );
    C(3*(i-1)+1:i*3,:) = simplify(res);
end
C = simplify(C);
C1 = C(1:3,:)
C2 = C(4:6,:)
C3 = C(7:9,:)

dq = [dq1,dq2,dq3];

% 1st factorization:

S = simplify([dq*C1; dq*C2; dq*C3])

skew_S = simplify(dM - 2*S)

% Skew-symmetry check

-skew_S == skew_S.'

% check energy conservation

simplify(dq*skew_S*dq.')

% 2nd factorization

S_prime = S - [0 -dq3 dq2; -dq3 0 dq1; dq2 -dq1 0] %[0 0 0; dq2 -dq1 0; -dq3 0 dq1]

skew_Sprime = simplify(dM - 2*S_prime)

% Skew-symmetry check

-skew_Sprime == skew_Sprime.'

% check energy conservation

simplify(dq*skew_Sprime*dq.')

% 3rd factorization

S_dprime = S - [0 0 0; -dq3 0 dq1; dq2 -dq1 0] %[0 0 0; dq2 -dq1 0; -dq3 0 dq1]

skew_Sdprime = simplify(dM - 2*S_dprime)

% Skew-symmetry check

-skew_Sdprime == skew_Sdprime.'

% check energy conservation

simplify(dq*skew_Sdprime*dq.')

% Linear parameterization :

tau = simplify(M*[ddq1; ddq2; ddq3] + c)

Y = simplify(jacobian(tau, [a1,a2,a3,a4,a5,a6]))

%% Functions:

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