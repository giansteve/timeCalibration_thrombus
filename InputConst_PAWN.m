function [ED,INPUT] = InputConst_PAWN(k_NonConst,Ns_const,konst_C,max_konstC)
%MORPH_INPUTCONST Create input structure in UQlab
% This function is meant to create the UQlab input structure that is then
% needed to create the Polynomial Chaos Expansion and perform the
% Senstivity Analysis.
% This version is designed for the "TimeCalibration" project.

%% Input
if ~ismember(1,k_NonConst(:))
    input.Marginals(1).Name = 'Dc';
    input.Marginals(1).Type = 'Uniform';
    input.Marginals(1).Parameters = [1.00e-10 1.00e-6];
else
    input.Marginals(1).Name = 'Dc';
    input.Marginals(1).Type = 'Constant';
    input.Marginals(1).Parameters = min([1.00e-10 1.00e-6]) + konst_C*(max([1.00e-10 1.00e-6]) - min([1.00e-10 1.00e-6]))/max_konstC;
end
if ~ismember(2,k_NonConst(:))
    input.Marginals(2).Name = 'k_c';
    input.Marginals(2).Type = 'Uniform';
    input.Marginals(2).Parameters = [200e2 200e4];
else
    input.Marginals(2).Name = 'k_c';
    input.Marginals(2).Type = 'Constant';
    input.Marginals(2).Parameters = min([200e2 200e4]) + konst_C*(max([200e2 200e4]) - min([200e2 200e4]))/max_konstC;
end
if ~ismember(3,k_NonConst(:))
    input.Marginals(3).Name = 'k_BP';
    input.Marginals(3).Type = 'Uniform';
    input.Marginals(3).Parameters = [8.000e-12 8.000e-8];
else
    input.Marginals(3).Name = 'k_BP';
    input.Marginals(3).Type = 'Constant';
    input.Marginals(3).Parameters = min([8.000e-12 8.000e-8]) + konst_C*(max([8.000e-12 8.000e-8]) - min([8.000e-12 8.000e-8]))/max_konstC;
end
if ~ismember(4,k_NonConst(:))
    input.Marginals(4).Name = 'C_t';
    input.Marginals(4).Type = 'Uniform';
    input.Marginals(4).Parameters = [10.00e2 10.00e4];
else
    input.Marginals(4).Name = 'C_t';
    input.Marginals(4).Type = 'Constant';
    input.Marginals(4).Parameters = min([10.00e2 10.00e4]) + konst_C*(max([10.00e2 10.00e4]) - min([10.00e2 10.00e4]))/max_konstC;
end
if ~ismember(5,k_NonConst(:))
    input.Marginals(5).Name = 'c_BPt';
    input.Marginals(5).Type = 'Uniform';
    input.Marginals(5).Parameters = [10 250e4];
else
    input.Marginals(5).Name = 'c_BPt';
    input.Marginals(5).Type = 'Constant';
    input.Marginals(5).Parameters = min([10 250e4]) + konst_C*(max([10 250e4]) - min([10 250e4]))/max_konstC;
end
if ~ismember(6,k_NonConst(:))
    input.Marginals(6).Name = 'T_Rt';
    input.Marginals(6).Type = 'Uniform';
    input.Marginals(6).Parameters = [1.00e-01 3.00e+00];
else
    input.Marginals(6).Name = 'T_Rt';
    input.Marginals(6).Type = 'Constant';
    input.Marginals(6).Parameters = min([1.00e-01 3.00e+00]) + konst_C*(max([1.00e-01 3.00e+00]) - min([1.00e-01 3.00e+00]))/max_konstC;
end
if ~ismember(7,k_NonConst(:))
    input.Marginals(7).Name = 'k_c,wall';
    input.Marginals(7).Type = 'Uniform';
    input.Marginals(7).Parameters = [100 10.00e4];
else
    input.Marginals(7).Name = 'k_c,wall';
    input.Marginals(7).Type = 'Constant';
    input.Marginals(7).Parameters = min([100 10.00e4]) + konst_C*(max([100 10.00e4]) - min([100 10.00e4]))/max_konstC;
end
if ~ismember(8,k_NonConst(:))
    input.Marginals(8).Name = 'c_AP';
    input.Marginals(8).Type = 'Uniform';
    input.Marginals(8).Parameters = [1.5e14/20 4.5e14/20];
else
    input.Marginals(8).Name = 'c_AP';
    input.Marginals(8).Type = 'Constant';
    input.Marginals(8).Parameters = min([1.5e14/20 4.5e14/20]) + konst_C*(max([1.5e14/20 4.5e14/20]) - min([1.5e14/20 4.5e14/20]))/max_konstC;
end
if ~ismember(9,k_NonConst(:))
    input.Marginals(9).Name = 'c_BPbt';
    input.Marginals(9).Type = 'Uniform';
    input.Marginals(9).Parameters = [100 250e3];
else
    input.Marginals(9).Name = 'c_BPbt';
    input.Marginals(9).Type = 'Constant';
    input.Marginals(9).Parameters = min([100 250e3]) + konst_C*(max([100 250e3]) - min([100 250e3]))/max_konstC;
end

if ~ismember(10,k_NonConst(:))
    input.Marginals(10).Name = 'Gammadot';
    input.Marginals(10).Type = 'Uniform';
    input.Marginals(10).Parameters = [0.1 50];
else
    input.Marginals(10).Name = 'Gammadot';
    input.Marginals(10).Type = 'Constant';
    input.Marginals(10).Parameters = min([0.1 50]) + konst_C*(max([0.1 50]) - min([0.1 50]))/max_konstC;
end
% M = size(input.Marginals,2);
try
    INPUT = uq_createInput(input);
catch
    uqlab
    INPUT = uq_createInput(input);
end
ED = uq_getSample(Ns_const,'lhs');
end % function