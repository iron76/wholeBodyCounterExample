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

S = solve(Nw_1x, Nw_1y, Eu_1, ...
    Nw_2x, Nw_2y, Eu_2, ...
    Nw_3x, Nw_3y, Eu_3, ...
    Nw_4x, Nw_4y, Eu_4, ...
    'a1x==0', 'a1y==0', 'o1==0', ...
    'a2x==0', 'a2y==0', 'o2==0', ...
    'a3x==0', 'a3y==0', 'o3==0', ...
    'a4x==0', 'a4y==0', 'o4==0', ...
    'tau1 == tau', 'tau2 == 0', 'tau3 == -tau',...
    f2x, f2y, f3x, f3y, f1x, f1y, ...
    fRx, fRy, fLx, fLy, tauR, tauL, ...
    a2x, a2y, o2, a3x, a3y, o3, ...
    a1x, a1y, o1, a4x, a4y, o4, ...
    tau1, tau2, tau3);
    
