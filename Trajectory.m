function [BStar,Increm] = Trajectory(N,M,GridLevel,PDF)
% This function computes the trajectories in the input space. This
% procedure is then used in the calculation o the Elementary Effects
% (Morris theory), as a preliminary Sensitivity Analysis. It is
% particularly suitable for model with a high multi-dimensionality.
%
% INPUT
% - N: number of trajectories to compute
% - M: dimension input space
% - GridLevel: numb. of subdivision into each single input dimension
%
% OUTPUT
% - BStar: matrix containing the evaluation points. BStar has dimensions
%           (M+1 x M x N)
% - Increm: incremental distance between points in input space
%
% Reference:
% 1. Saltelli et al. (2008), Global Sensitivity Analisys: the Primer

message = 'Computing trajectories ... \n';
fprintf(message)
time = tic;

%% 2.1 - Computing deltas
Delta = GridLevel / (2*GridLevel - 2);
BStar = zeros(M+1,M,N);

for BS = 1:N
    
    for InpIdx = 1:M
        DeltaToEE(InpIdx) = Delta;
    end
    
    %% 2.2 - Building Experimantal Design
    % Assuming PDF for the input dimensions
    PDFInput = PDF.PDFInput;
    mu = PDF.mu;
    sigma = PDF.sigma;
    truncations = PDF.truncations;
    UnifBounds = zeros(M,2);    
    testNorm = 'normal';
    InputGrid = zeros(M,GridLevel);
    
    for InpIdx = 1:M
        
        if PDFInput(InpIdx) == testNorm % normal PDF
            PDFs(InpIdx) = truncate(makedist(testNorm,mu(InpIdx),...
                sigma(InpIdx)),truncations(InpIdx,1),truncations(InpIdx,2));
            DeltaToEE(InpIdx) = 1/(GridLevel-1);
            quantiles = linspace(0,1,GridLevel);
            InputGrid(InpIdx,:) = icdf(PDFs(InpIdx),quantiles);
        elseif truncations(InpIdx,1) == 0 && truncations(InpIdx,2) == 1 % uniform PDF [0,1]
            UnifBounds(InpIdx,1:2) = [0,1];
            InputGrid(InpIdx,:) = linspace(0,1,GridLevel);
        else % uniform PDF ~= [0,1]
            UnifBounds(InpIdx,1:2) = [truncations(InpIdx,1) ...
                truncations(InpIdx,2)];
            DeltaToEE(InpIdx) = DeltaToEE(InpIdx)*(truncations(InpIdx,2)-...
                truncations(InpIdx,1));
            InputGrid(InpIdx,:) = linspace(truncations(InpIdx,1)...
                ,truncations(InpIdx,2),GridLevel);
        end
    end
    
    %% 2.3 - x star position vector
    % x star is randomly chosen base value of the Experimental Design
    
    % TEST! to avoid points outside the input domain
    xStar = zeros(M,1);
    for XIdx = 1:M
        if PDFInput(XIdx) == testNorm
            testMin = icdf(PDFs(XIdx),(quantiles) + (1/(GridLevel-1)));
            testLogic = testMin < max(InputGrid(XIdx,:));
            AmountAvailNorm = randi(sum(testLogic));
            xStar(XIdx,:) = InputGrid(XIdx,AmountAvailNorm);
        else
            test = InputGrid(XIdx,:) + DeltaToEE(XIdx);
            testLogic = test < truncations(XIdx,2);
            AmountAvailUnif = randi(sum(testLogic));
            xStar(XIdx,:) = InputGrid(XIdx,AmountAvailUnif);
        end
    end
    xStar = xStar';
    
    %% 2.4 BStar matrices
    % D star matrix
    DStar = zeros(M,M);
    DSinput = [1 -1];
    for DSIdx = 1:M
        DStar(DSIdx,DSIdx) = DSinput(randi(2));
    end
    
    % B matrix
    BTemp = ones(M+1,M);
    BTemp = tril(BTemp);
    B = [zeros(1,M) ; BTemp];
    B = B(1:M+1,:);
    
    % J matrix
    JMatrix = ones(M+1,M);
    
    % P star matrix
    % possibility to randomize the order of points for the evaluation.
    % Suggested only in case of i.i.d. random variables
    PStar = eye(M);
    
    % Construction of BStar matrices
    Increm = ((1/2)*((2*B-JMatrix)*DStar + JMatrix));
    for InpIdx = 1:M
        if PDFInput(InpIdx) == testNorm
            QuantFin = round(cdf(PDFs(InpIdx),xStar(InpIdx)),4) + DeltaToEE(InpIdx);
            ValFromQuant = icdf(PDFs(InpIdx),QuantFin);
            Increm(:,InpIdx) = Increm(:,InpIdx) * ValFromQuant;
        else
            Increm(:,InpIdx) = Increm(:,InpIdx) * DeltaToEE(InpIdx);
        end
    end
    
    for InpIdx = 1:M
        BStar(:,:,BS) = Increm+xStar;
    end
    
    BStar(:,:,BS) = BStar(:,:,BS) * PStar;
    
end

%%
% Print out on Command Window the total amount of time used for the
% computation of Trajectories
time = toc(time);
message = 'Time computing trajectories: %4.2f seconds \n';
fprintf(message,time)

%%
% Elementary Effect Evaluation interfaced with COMSOL software
%
% Author: Gian Marco Melito         Date: December 2018
% Project: Aorta.COMSOL
% TU Graz
