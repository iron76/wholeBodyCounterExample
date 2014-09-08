% Genova, 13 July 2014
%
% We consider the simple counter example shown in the Frontiers paper. The
% example shows that the GCoP within the support polygon is not a
% sufficient condition for the local FRI to be stable. The geometry of the
% system is the one represented in the figure zmpCounterExample.pdf.

clear all
close all
clc

syms m1 m2 m3 m4 tau1 tau2 tau3 tau h l real
syms fRx fLx tauR tauL f1x f2x f3x real  
syms fRy fLy tauR tauL f1y f2y f3y real 
syms a1x a1y o1 m1 real
syms a2x a2y o2 m2 real
syms a3x a3y o3 m3 real
syms a4x a4y o4 m4 real
syms I1 I2 I3 I4 real

% Newton on left foot
% m1 a1x = fLx - f1x
% m1 a1y = fLy - f1y - m1 * g
% Euler  on right link
% I1 o1  =  tauL - tau1 + f1x*h/2 + fLx*h/2

Nw_1x = 'fLx - f1x == m1 * a1x';
Nw_1y = 'fLy - f1y - m1 * g == m1 * a1y';
Eu_1  = ' tauL - tau1 + f1x*h/2 + fLx*h/2 == I1 * o1';

% Newton on left link
% m2 a2x = f1x + f2x
% m2 a2y = f1y + f2y - m2 * g
% Euler  on right link
% I2 o2  = f1x * l + f2y * l + tau1 + tau2

Nw_2x = 'f1x + f2x == m2 * a2x';
Nw_2y = 'f1y + f2y - m2 * g == m2 * a2y';
Eu_2  = 'f1x * l + f2y * l + tau1 + tau2 == I2 * o2';

% Newton on right link
% m3 a3x = f3x - f2x 
% m3 a3y = f3y - f2y - m3g
% Euler  on right link
% I3 o3  = f3x * l + f2y * l - tau2 + tau3

Nw_3x = 'f3x - f2x == m3 * a3x';
Nw_3y = 'f3y - f2y - m3 * g == m3 * a3y';
Eu_3  = 'f3x * l + f2y * l - tau2 + tau3 == I3 * o3';

% Newton on right foot
% m4 a4x = fRx - f3x
% m4 a4y = fRy - f3y - m4 * g
% Euler  on right link
% I4 o4  =  tauR - tau3 + f3x*h/2 + fRx*h/2

Nw_4x = 'fRx - f3x == m4 * a4x';
Nw_4y = 'fRy - f3y - m4 * g == m4 * a4y';
Eu_4  = ' tauR - tau3 + f3x*h/2 + fRx*h/2 == I4 * o4';


% Kinematic constraints on the foot
% rotation. Feet are rotating around
% the edge. The only free variables are
% d2qR and d2qL, corresponding to the left
% and right foot tipping acceleration.

syms    qL    q1    q2    q3 real
syms dth_L dth_1 dth_2 dth_3 real

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

qL    = 3*pi/4;     q1 = -pi/2;     q2 = -pi/2;     q3 = -pi/2;
dth_L = 3*pi/4;  dth_1 =  pi/4;  dth_2 =  pi/4;  dth_3 =  pi/4;

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

% Angular velocities
Kn_1o = 'o1 == d2qL';
Kn_2o = 'o2 == d2qL + d2q1';
Kn_3o = 'o3 == d2qL + d2q1 + d2q2';
Kn_4o = 'o4 == d2qL + d2q1 + d2q2 + d2q3';

%Torques at the feet
Ct_Lf =  'tauL ==  h*fRy';
Ct_Rf =  'tauR == -h*fLy';

S = solve(...
Nw_1x, Nw_1y, Eu_1, ...   3eq
Nw_2x, Nw_2y, Eu_2,  ...  6eq
Nw_3x, Nw_3y, Eu_3,  ...  9eq
Nw_4x, Nw_4y, Eu_4,  ... 12eq
Kn_1x, Kn_1y, Kn_1o, ... 15eq
Kn_2x, Kn_2y, Kn_2o, ... 18eq
Kn_3x, Kn_3y, Kn_3o, ... 21eq
Kn_4x, Kn_4y, Kn_4o, ... 24eq
Kn_Rx, Kn_Ry,        ... 26eq
Ct_Lf, Ct_Rf,        ... 28eq
'tau1 == tau', 'tau2 == 0', 'tau3 == -tau',... 31eq
d2qL, d2q1, d2q2, d2q3, ...         4kn
f2x, f2y, f3x, f3y, f1x, f1y, ...   10kn
fRx, fRy, fLx, fLy, tauR, tauL, ... 16kn
a2x, a2y, o2, a3x, a3y, o3, ...     22kn
a1x, a1y, o1, a4x, a4y, o4, ...     28kn
tau1, tau2, tau3);                 %31kn
    
