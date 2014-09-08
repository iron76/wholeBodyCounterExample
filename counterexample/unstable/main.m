% Genova, 13 July 2014
%
% We consider the simple counter example shown in the Frontiers paper. The
% example shows that the GCoP within the support polygon is not a
% sufficient condition for the local FRI to be stable. The geometry of the
% system is the one represented in the figure zmpCounterExample.pdf.

clear all
close all
clc

syms g h l dx positive
syms tau1 tau2 tau3 tau real
syms fRx fLx tauR tauL f1x f2x f3x real  
syms fRy fLy tauR tauL f1y f2y f3y real 
syms a1x a1y o1  real
syms a2x a2y o2  real
syms a3x a3y o3  real
syms a4x a4y o4  real
syms I1 I2 I3 I4 I positive
syms m1 m2 m3 m4 m positive

% Kinematic constraints on the foot
% rotation. Feet are rotating around
% the edge. The only free variables are
% d2qR and d2qL, corresponding to the left
% and right foot tipping acceleration.

syms    qL q1  q2  q3 real
syms dth_L dth_1 dth_2 dth_3 real

% qR = qL + q1 + q2 + q3 + pi

xcom1 = sqrt(2)*h*cos(qL) + h/2 * cos(qL+dth_L);
ycom1 = sqrt(2)*h*sin(qL) + h/2 * sin(qL+dth_L);
xcom2 = sqrt(2)*h*cos(qL) + l*cos(qL+q1+dth_1);
ycom2 = sqrt(2)*h*sin(qL) + l*sin(qL+q1+dth_1);
xcom3 = sqrt(2)*h*cos(qL) + sqrt(2)*l*cos(qL+q1) + l*cos(qL+q1+q2+dth_2);
ycom3 = sqrt(2)*h*sin(qL) + sqrt(2)*l*sin(qL+q1) + l*sin(qL+q1+q2+dth_2);
xcom4 = sqrt(2)*h*cos(qL) + sqrt(2)*l*cos(qL+q1) + sqrt(2)*l*cos(qL+q1+q2) + h/2*cos(qL+q1+q2+q3+dth_3);
ycom4 = sqrt(2)*h*sin(qL) + sqrt(2)*l*sin(qL+q1) + sqrt(2)*l*sin(qL+q1+q2) + h/2*sin(qL+q1+q2+q3+dth_3);
xR    = sqrt(2)*h*cos(qL) + sqrt(2)*l*cos(qL+q1) + sqrt(2)*l*cos(qL+q1+q2) + sqrt(2)*h*cos(qL+q1+q2+q3); 
yR    = sqrt(2)*h*sin(qL) + sqrt(2)*l*sin(qL+q1) + sqrt(2)*l*sin(qL+q1+q2) + sqrt(2)*h*sin(qL+q1+q2+q3);


Jx1    = jacobian(xcom1, [qL    q1    q2    q3]);
Jy1    = jacobian(ycom1, [qL    q1    q2    q3]);
Jx2    = jacobian(xcom2, [qL    q1    q2    q3]);
Jy2    = jacobian(ycom2, [qL    q1    q2    q3]);
Jx3    = jacobian(xcom3, [qL    q1    q2    q3]);
Jy3    = jacobian(ycom3, [qL    q1    q2    q3]);
Jx4    = jacobian(xcom4, [qL    q1    q2    q3]);
Jy4    = jacobian(ycom4, [qL    q1    q2    q3]);
JxR    = jacobian(   xR, [qL    q1    q2    q3]);
JyR    = jacobian(   yR, [qL    q1    q2    q3]);

dth_L = 3*pi/4;  dth_1 =  pi/4;  dth_2 =  pi/4;   dth_3 =  pi/4;    dth_R = 5*pi/4;
qL    = 3*pi/4;     q1 = -pi/2;     q2 = -pi/2;     q3  = -pi/2;       qR =  pi/4;
%

Jx1n    = vpa(eval(Jx1), 4);
Jy1n    = vpa(eval(Jy1), 4);
Jx2n    = vpa(eval(Jx2), 4);
Jy2n    = vpa(eval(Jy2), 4);
Jx3n    = vpa(eval(Jx3), 4);
Jy3n    = vpa(eval(Jy3), 4);
Jx4n    = vpa(eval(Jx4), 4);
Jy4n    = vpa(eval(Jy4), 4);
JxRn    = vpa(eval(JxR), 4);
JyRn    = vpa(eval(JyR), 4);

Jy2n = [ - 1.0*h, 0, 0, 0];
Jy4n = [ 2.0*l - 1.0*h, 2.0*l, l, 0];

% joint accelerations
syms d2qL d2q1 d2q2 d2q3 real

% kinematic constraints
a1x_q = Jx1n * [d2qL d2q1 d2q2 d2q3]';
a1y_q = Jy1n * [d2qL d2q1 d2q2 d2q3]';
a2x_q = Jx2n * [d2qL d2q1 d2q2 d2q3]';
a2y_q = Jy2n * [d2qL d2q1 d2q2 d2q3]';
a3x_q = Jx3n * [d2qL d2q1 d2q2 d2q3]';
a3y_q = Jy3n * [d2qL d2q1 d2q2 d2q3]';
a4x_q = Jx4n * [d2qL d2q1 d2q2 d2q3]';
a4y_q = Jy4n * [d2qL d2q1 d2q2 d2q3]';
aRx_q = JxRn * [d2qL d2q1 d2q2 d2q3]';
aRy_q = JyRn * [d2qL d2q1 d2q2 d2q3]';

Kn_1x = ['a1x == ' char(a1x_q)];
Kn_1y = ['a1y == ' char(a1y_q)];
Kn_2x = ['a2x == ' char(a2x_q)];
Kn_2y = ['a2y == ' char(a2y_q)];
Kn_3x = ['a3x == ' char(a3x_q)];
Kn_3y = ['a3y == ' char(a3y_q)];
Kn_4x = ['a4x == ' char(a4x_q)];
Kn_4y = ['a4y == ' char(a4y_q)];
Kn_Rx = ['0 == '   char(aRx_q)];
Kn_Ry = ['0 == '   char(aRy_q)];

% d2q2 = d2qR - d2qL - d2q1 - d2q3;

% Angular velocities
o1_q = d2qL;
o2_q = d2qL + d2q1;
o3_q = d2qL + d2q1 + d2q2;
o4_q = d2qL + d2q1 + d2q2 + d2q3;

%Torques at the feet
Ct_Lf =  'tauL - h*fRy == 0';
Ct_Rf =  'tauR + h*fLy == 0';

% Newton on left foot
% m1 a1x = fLx - f1x
% m1 a1y = fLy - f1y - m1 * g
% Euler  on right link
% I1 o1  =  tauL - tau1 + f1x*h/2 + fLx*h/2

Nw_1x = ['fLx - f1x ==' char(m1 * a1x_q)];
Nw_1y = ['fLy - f1y - m1 * g ==' char(m1 * a1y_q)];
Eu_1  = [' tauL - tau1 + f1x*h/2 + fLx*h/2 ==' char(I1 * o1_q)];

% Newton on left link
% m2 a2x = f1x + f2x
% m2 a2y = f1y + f2y - m2 * g
% Euler  on right link
% I2 o2  = f1x * l + f2y * l + tau1 + tau2

Nw_2x = ['f1x + f2x ==' char(m2 * a2x_q)];
Nw_2y = ['f1y + f2y - m2 * g ==' char(m2 * a2y_q)];
Eu_2  = ['f1x * l + f2y * l + tau1 + tau2 == ' char(I2 * o2_q)];

% Newton on right link
% m3 a3x = f3x - f2x 
% m3 a3y = f3y - f2y - m3g
% Euler  on right link
% I3 o3  = f3x * l + f2y * l - tau2 + tau3

Nw_3x = ['f3x - f2x == ' char(m3 * a3x_q)];
Nw_3y = ['f3y - f2y - m3 * g == ' char(m3 * a3y_q)];
Eu_3  = ['f3x * l + f2y * l - tau2 + tau3 == ' char(I3 * o3_q)];

% Newton on right foot
% m4 a4x = fRx - f3x
% m4 a4y = fRy - f3y - m4 * g
% Euler  on right link
% I4 o4  =  tauR - tau3 + f3x*h/2 + fRx*h/2

Nw_4x = ['fRx - f3x == ' char(m4 * a4x_q)];
Nw_4y = ['fRy - f3y - m4 * g == ' char(m4 * a4y_q)];
Eu_4  = ['tauR - tau3 + f3x*h/2 + fRx*h/2 == ' char(I4 * o4_q)];

% syms _d2qL_ _d2q1_ _d2q2_ _d2q3_ real

S = solve(...
Nw_1x, Nw_1y, Eu_1, ...   3eq
Nw_2x, Nw_2y, Eu_2,  ...  6eq
Nw_3x, Nw_3y, Eu_3,  ...  9eq
Nw_4x, Nw_4y, Eu_4,  ... 12eq
Kn_Rx, Kn_Ry,        ... 14eq
Ct_Lf, Ct_Rf,        ... 16eq
'tau1 == tau', 'tau2 == 0', 'tau3 == -tau',... 19eq
d2qL, d2q1, d2q2, d2q3, ...         4kn
f2x, f2y, f3x, f3y, f1x, f1y, ...   10kn
fRx, fRy, fLx, fLy, tauR, tauL, ... 16kn
tau1, tau2, tau3, 'IgnoreAnalyticConstraints', false, 'Real', true, 'PrincipalValue', true);                % 19kn
    
% S0 = solve(...
% Nw_1x, Nw_1y, Eu_1, ...   3eq
% Nw_2x, Nw_2y, Eu_2,  ...  6eq
% Nw_3x, Nw_3y, Eu_3,  ...  9eq
% Nw_4x, Nw_4y, Eu_4,  ... 12eq
% 'tau1 == tau', 'tau2 == 0', 'tau3 == -tau',... 15eq
% 'd2qL == 0'  , 'd2q1 == 0', 'd2q2 == 0'   , 'd2q3 ==0', ...         19eq
% d2qL, d2q1, d2q2, d2q3, ...         4kn
% f2x, f2y, f3x, f3y, f1x, f1y, ...   10kn
% fRx, fRy, fLx, fLy, tauR, tauL, ... 16kn
% tau1, tau2, tau3);                % 19kn
% 

d2qL = S.d2qL;
d2q1 = S.d2q1;
d2q2 = S.d2q2;
d2q3 = S.d2q3;
d2qR = d2q2 + d2qL + d2q1 + d2q3;


fLx  = S.fLx;
fLy  = S.fLy;
f1x  = S.f1x;
f1y  = S.f1y;
f2x  = S.f2x;
f2y  = S.f2y;
f3x  = S.f3x;
f3y  = S.f3y;
fRx  = S.fRx;
fRy  = S.fRy;

tauR = S.tauR;
tauL = S.tauL;

tau1 = S.tau1;
tau2 = S.tau2;
tau3 = S.tau3;

% simple(eval(fLx - f1x - m1 * a1x_q))
% simple(eval(f3x * l + f2y * l - tau2 + tau3 - I3 * o3_q))
% simple(eval(d2q1 + d2q3))


dqL_val =  1;
dqR_val =  -1;
% dqR = dqL + dq1 + dq2 + dq3
% 0   = dqL + dq1 + dq2 + dq3 - dqR
J = [JxRn 0; JyRn 0; 1 1 1 1 -1];
dq_val = -J(:, [2 3 4])^(-1)*J(:, [1 5])*[dqL_val, dqR_val]';
simple(dq_val(1) - dq_val(3))

% Jr = [JxRn; JyRn];
% dq_val = -Jr(:, [2 4])^(-1)*Jr(:, [1 3])*[dqL_val, 0]';

% I  = m*l^2;
I1 = I;   I2 = I;   I3 = I;   I4 = I;
m1 = m/4; m2 = m/4; m3 = m/4; m4 = m/4;
h  = l/2;
%tau = 1/6*g*l*m;

FRI_r_eq = tauR/fRy;
FRI_l_eq = tauL/fLy;

x1       = sqrt(2)*h*cos(qL);
y1       = sqrt(2)*h*sin(qL);
FRI_l    = (tau1-y1*f1x)/(f1y+m*g/4);

x3       = sqrt(2)*h*cos(qL) + sqrt(2)*l*cos(qL+q1) + sqrt(2)*l*cos(qL+q1+q2);
y3       = sqrt(2)*h*sin(qL) + sqrt(2)*l*sin(qL+q1) + sqrt(2)*l*sin(qL+q1+q2);
FRI_r    = (tau3-y3*f3x)/(f3y+m*g/4);


simple(eval(d2qR+d2qL))
simple(eval(fLx+fRx))
simple(eval(fLy-fRy))
simple(eval(FRI_l-l/2))
simple(eval(FRI_r+l/2))
