function Y = timeCal_RAmodel(X,P)
% model for the evaluation of the PCE for the Reliability Analysis problem


% evaluate PCE
output_PCE = uq_evalModel(P.mod,X);

% correct output by threshold
Y = zeros(size(X,1),1) - 1;
Y(abs(output_PCE) <= P.threshold) = 1;

end % function