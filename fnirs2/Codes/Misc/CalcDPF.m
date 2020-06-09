 
function [DPF] = CalcDPF (A,lambda)
 
% Parameters
 
a =  223.3;
b =  0.05624;
c = 0.8493;
d = -5.723E-7;
e =  0.001245;
f = -0.9025;
 
DPF = a + b.*A.^c + d.*lambda.^3 + e.*lambda.^2 + f.*lambda;
%DPF = sprintf('%.2f\n',DPF);