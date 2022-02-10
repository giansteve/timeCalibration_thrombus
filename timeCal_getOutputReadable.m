function [outToPCE] = timeCal_getOutputReadable(OUTPUT)
%MORPH_GETOUTPUT_READABLE This function transforms the outputs computed
%from OpenFOAM into something that is readable for the computation of the
%PCE
% INPUT:
%       - OUTPUT: the structure coming out from the reading of the
%           simulation results

for nsim = 1:length(OUTPUT)
    
    %% Time
    outToPCE.time(:,nsim) = OUTPUT(nsim).time(20:end);
    
    %% H/S
    %     logic_temp = OUTPUT(nsim).H_S.phi_c(:,2:end) >= thresholdThr;
    %     OUTPUT(nsim).H_S.output = (sum(logic_temp,2) / OUTPUT(nsim).H_S.Tot_pts * OUTPUT(nsim).H_S.Length(2)) / OUTPUT(nsim).H_S.Length(2);
    outToPCE.H_S(:,nsim) = OUTPUT(nsim).H_S.output(20:end);
    
    %% L/S
    %     logic_temp = OUTPUT(nsim).L_S.phi_c(:,2:end) >= thresholdThr;
    %     OUTPUT(nsim).L_S.output = (sum(logic_temp,2) / OUTPUT(nsim).L_S.Tot_pts * OUTPUT(nsim).L_S.Length(1)) / OUTPUT(nsim).H_S.Length(2);
    outToPCE.L_S(:,nsim) = OUTPUT(nsim).L_S.output(20:end);
    
    %% Surface area
    %     outputToPCE.SA(:,n_ns) = OUTPUT(n_ns).SA.output(1000:50:end);
    
    %% Volume
    %     outputToPCE.Vol(:,n_ns) = OUTPUT(n_ns).Vol.output(1000:50:end).*1.075e-06;
    
end


end

