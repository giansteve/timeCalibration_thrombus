function [INPUT,exp_design] = timeCal_createInput(Ns,nonConst,expDesign_switch,ED,Eplot)
%CREATEINPUT_MORPHOLOGICAL Create input structure in UQlab
% This function is meant to create the UQlab input structure that is then
% needed to create the Polynomial Chaos Expansion and perform the
% Senstivity Analysis.
% This version is designed for the "TimeCalibration" project.

root_destination = pwd;

if isempty(ED)
    kbp = [8.000e-12 8.000e-8];
    cbpt = [10 250e4];
    
    %% Input
    if ismember(1,nonConst(:))
        input.Marginals(1).Name = '$D_{c}$';
        input.Marginals(1).Type = 'Uniform';
        input.Marginals(1).Parameters = [1.00e-10 1.00e-6];
    else
        input.Marginals(1).Name = '$D_{c}$';
        input.Marginals(1).Type = 'Constant';
        input.Marginals(1).Parameters = mean([1.00e-10 1.00e-6]);
    end
    
    if ismember(2,nonConst(:))
        input.Marginals(2).Name = '$k_{c}$';
        input.Marginals(2).Type = 'Uniform';
        input.Marginals(2).Parameters = [200e2 200e4];
    else
        input.Marginals(2).Name = '$k_{c}$';
        input.Marginals(2).Type = 'Constant';
        input.Marginals(2).Parameters = mean([200e2 200e4]);
    end
    
    if ismember(3,nonConst(:))
        input.Marginals(3).Name = '$k_{BP}$';
        input.Marginals(3).Type = 'Uniform';
        input.Marginals(3).Parameters = kbp;
    else
        input.Marginals(3).Name = '$k_{BP}$';
        input.Marginals(3).Type = 'Constant';
        input.Marginals(3).Parameters = mean(kbp);
    end
    
    if ismember(4,nonConst(:))
        input.Marginals(4).Name = '$C_{t}$';
        input.Marginals(4).Type = 'Uniform';
        input.Marginals(4).Parameters = [10.00e2 10.00e4];
    else
        input.Marginals(4).Name = '$C_{t}$';
        input.Marginals(4).Type = 'Constant';
        input.Marginals(4).Parameters = mean([10.00e2 10.00e4]);
    end
    
    if ismember(5,nonConst(:))
        input.Marginals(5).Name = '$c_{BPt}$';
        input.Marginals(5).Type = 'Uniform';
        input.Marginals(5).Parameters = cbpt;
    else
        input.Marginals(5).Name = '$c_{BPt}$';
        input.Marginals(5).Type = 'Constant';
        input.Marginals(5).Parameters = mean(cbpt);
    end
    
    if ismember(6,nonConst(:))
        input.Marginals(6).Name = '$T_{Rt}$';
        input.Marginals(6).Type = 'Uniform';
        input.Marginals(6).Parameters = [1.00e-01 3.00e+00];
    else
        input.Marginals(6).Name = '$T_{Rt}$';
        input.Marginals(6).Type = 'Constant';
        input.Marginals(6).Parameters = mean([1.00e-01 3.00e+00]);
    end
    
    if ismember(7,nonConst(:))
        input.Marginals(7).Name = '$k_{c,wall}$';
        input.Marginals(7).Type = 'Uniform';
        input.Marginals(7).Parameters = [100 10.00e4];
    else
        input.Marginals(7).Name = '$k_{c,wall}$';
        input.Marginals(7).Type = 'Constant';
        input.Marginals(7).Parameters = mean([100 10.00e4]);
    end
    
    if ismember(8,nonConst(:))
        input.Marginals(8).Name = '$c_{AP}$';
        input.Marginals(8).Type = 'Uniform';
        input.Marginals(8).Parameters = [1.5e14/20 4.5e14/20];
    else
        input.Marginals(8).Name = '$c_{AP}$';
        input.Marginals(8).Type = 'Constant';
        input.Marginals(8).Parameters = mean([1.5e14/20 4.5e14/20]);
    end
    
    if ismember(9,nonConst(:))
        input.Marginals(9).Name = '$c_{BPbt}$';
        input.Marginals(9).Type = 'Uniform';
        input.Marginals(9).Parameters = [100 250e3];
    else
        input.Marginals(9).Name = '$c_{BPbt}$';
        input.Marginals(9).Type = 'Constant';
        input.Marginals(9).Parameters = mean([100 250e3]);
    end
    
    if ismember(10,nonConst(:))
        input.Marginals(10).Name = '$gamma_{t}$';
        input.Marginals(10).Type = 'Uniform';
        input.Marginals(10).Parameters = [0.1 50];
    else
        input.Marginals(10).Name = '$gamma_{t}$';
        input.Marginals(10).Type = 'Constant';
        input.Marginals(10).Parameters = mean([0.1 50]);
    end
    
    M = size(input.Marginals,2);
    try
        INPUT = uq_createInput(input);
    catch
        uqlab
        INPUT = uq_createInput(input);
    end
    
    if expDesign_switch == 1
        exp_design = uq_getSample(Ns,'lhs');
    end
    
else
    try
        iOpts.Inference.Data = ED;
        iOpts.Copula.Type = 'Independent';
        INPUT = uq_createInput(iOpts);
    catch
        uqlab
        iOpts.Inference.Data = ED;
        iOpts.Copula.Type = 'Independent';
        INPUT = uq_createInput(iOpts);
    end
end

if Eplot == 1
% make dist for output
try
    dest_plot = sprintf('Plot\1500sims');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end

% plot input
figure('Visible','off')
plotmatrix(ED)
GM_printBMP(500,500,sprintf('INPUTstats'))
GM_printEPS(500,500,sprintf('INPUTstats'))
close

cd(root_destination)
end





end

