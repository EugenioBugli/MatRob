clear
clc

syms q1 q2 q3 q4 real
syms a1 a2 a3 a4 a5 real
syms dc1 dc2 dc3 dc4 real
syms l1 l2 l3 l4 l real 
syms dq1 dq2 dq3 dq4 real
syms Ic1 Ic2 Ic3 Ic4 real
syms m1 m2 m3 m4 real
syms rc1x rc1y rc2x rc2y real
syms g0 real


% 2R Robot

M =simplify([Ic1 + Ic2 + l1^2*m1 + l2^2*m2 + m1*rc1x^2 + m1*rc1y^2 + m2*rc2x^2 + m2*rc2y^2 + 2*l1*m1*rc1x + 2*l2*m2*rc2x + l1^2*m2*cos(q2)^2 + l1^2*m2*sin(q2)^2 + 2*l1*l2*m2*cos(q2) + 2*l1*m2*rc2x*cos(q2) - 2*l1*m2*rc2y*sin(q2), m2*l2^2 + 2*m2*l2*rc2x + l1*m2*cos(q2)*l2 + m2*rc2x^2 + l1*m2*cos(q2)*rc2x + m2*rc2y^2 - l1*m2*sin(q2)*rc2y + Ic2;
                                                                                                     m2*l2^2 + 2*m2*l2*rc2x + l1*m2*cos(q2)*l2 + m2*rc2x^2 + l1*m2*cos(q2)*rc2x + m2*rc2y^2 - l1*m2*sin(q2)*rc2y + Ic2,                                                              m2*l2^2 + 2*m2*l2*rc2x + m2*rc2x^2 + m2*rc2y^2 + Ic2])

