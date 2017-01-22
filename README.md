# Singular-LQ
MATLAB Toolbox for Singular LQ Sliding Surface design.

% [T,SS,QC,Ueqn,algorithm] = SSLQ(A,B,Qb,Gain)
% SingularLQ: 1) Calculates transformation T from matrices A and B,  to a Controllable Canonical Form:
% z = Tx,
% 2) Depending on Weithing Matrix Qb, calculates the Sliding Surface and its
% r-1 time derivatives [s \dot{s} \ddot{s} ...]^{T} = SS*z, 
% 3) Generates Quasi-Continuous Algorithm of order r with the specified Gain,
% 4) Calculates nominal equivalent control ueqn = Ueqn*z, so the complete control law is u = ueqn
% + QC
