%% ELEMENTARY EFFECT EVALUATION (ElEv_Aorta.COMSOL)
% Evaluation o the Elementary Effect of a given function. This method is
% ideal when the number of input factor is too large to allow application
% of computationally expensive variance-based techniques.
%
% The aim is determining which input factor could be considered to have
% effects on the aoutput such as:
% - negligible
% - linear and additive
% - non linear or involved in interactions
%
% INPUT
% - M: dimension of the input space
% - N_seq: number of trajectories to be build into input space
% - GridLevel: numb. of subdivision into each single input dimension
%
% OUTPUT
% - Mu: overall influence of the variable on the output
% - MuStar: solves type II error; i.e. when the model in nonmonotonic or
%           has interactions
% - VarEE: estimates the ensemble of the variable'S effects whether
%          nonlinear and/or due to interactions
%
clc; clear;
timeT = tic;

set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\general')
addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\SA')
addpath('M:\IFM\User\melito\Server\Projects\matlab_funct')


%% 1 - Setting the evaluation

% N_seq = round(linspace(10,100,10));    % if need convergence trajectories

N_seq = 100;                             % number of trajectories
GridLevel = 50;                          % size grid in input axis

rvNames =           {'${D}_{\mathrm{c}}$','${k}_{\mathrm{c}}$','${k}_{\mathrm{{BP}}}$','${c}_{\mathrm{t}}$',...
                    '${c}_{\mathrm{{BPt}}}$','$\overline{T}_{\mathrm{Rt}}$','${k}_{\mathrm{{c,wall}}}$',...
                    '${c}_{\mathrm{{AP}}}$','${c}_{\mathrm{{BPbt}}}$','$\dot{\overline{\gamma}}_t$'};
PDF.PDFInput =         ['uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform';
                        'uniform'];
% insert if PDF is normal; otherwise 0
PDF.mu =            [0 0 0 ...
    0 0 0 ...
    0 0];
% insert if PDF is normal; otherwise 0
PDF.sigma =         [0 0 0 ...
    0 0 0 ...
    0 0];
% if uniform PDF write moments, otherwise insert truncations
PDF.truncations =   [1e-10,1e-06;
                    2.00e+04,2.00e+06;
                    8e-13 8e-8;
                    1.0e+03,1.0e+05;
                    1e+01,2.5e+06;
                    0.1 1;
                    100,10e4;
                    [1.5e14 4.5e14]./20;
                    100,2.5e+05;
                    0.1,20.0];
% PDF.truncations =   [0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1;
%     0 1];

%% 2 - Collecting Experimental Design
% computing the trajectories into the multi dimensional input space
M = size(rvNames,2);
for Nidx = 1:length(N_seq)
    N = N_seq(Nidx);
    numSim = N * (M+1);                 % tot numb simulations
    fprintf('   ---   ---   ---   ---   ---   ')
    SimuRun = '\nSimulation for %4d trajectories \n';
    fprintf(SimuRun,N)
    RunNeeded = 'Number of simulations ongoing: %4d \n';
    fprintf(RunNeeded,numSim)
    
    [BStar] = Trajectory(N,M,GridLevel,PDF); % Collecting trajectories
    
    %     %% 3 - Model definition and solution
    %     time = tic;
    %     source = pwd;
    %     Y = zeros(M+1,N);
    %     i = 0;
    %     for SolIdx = 1:N
    %         for SolIdx2 = 1:M+1
    %             i = i+1;
    %             Y(SolIdx2,SolIdx) = InputFileConstruction(BStar(SolIdx2,:,SolIdx),source);
    %         end
    %     end
    %     time = toc(time);
    %     message = 'Solver time: %4.2f seconds \n';
    %     fprintf(message,time)
    
    %     %% 4 - Evaluating EE indices
    %     [Mu,MuStar,VarEE] = EE_indices(BStar,Y,M,N);
    %     Solut(:,:,Nidx) = [Mu;MuStar;VarEE];
    %
end
fprintf('Writing in Excel ... \n')
excelFileName = 'EEeffect_simRun.xlsx';
for nseq = 1:N_seq
xlsappend(excelFileName,BStar(:,:,nseq),'sims')
end
% %% 5 - Plots
% if length(N_seq) == 1
%     BarPlot(Mu,MuStar,VarEE,rvNames)
% end
% %% 5* - Convergence assessment
% % In order to find out the minimum and feasible number of trajectories
% % needed for the Sensitivity Assessment, the study will repeat with
% % different number of trajectories
% if length(N_seq) > 1
%     for DimInp = 1:M
%         for jj = 1:length(N_seq)
%             MuPlot(jj,DimInp) = Solut(1,DimInp,jj);
%             MuStarPlot(jj,DimInp) = Solut(2,DimInp,jj);
%             VarEEPlot(jj,DimInp) = Solut(3,DimInp,jj);
%         end
%     end
%     
%     figure                      % convergence of MU
%     for DimInp = 1:M
%         semilogx(N_seq,MuPlot(:,DimInp),'-x')
%         hold on
%     end
%     hold off
%     title('\mu')
%     legend(rvNames)
%     
%     figure                      % convergence of MU*
%     for DimInp = 1:M
%         semilogx(N_seq,MuStarPlot(:,DimInp),'-x')
%         hold on
%     end
%     hold off
%     title('\mu*')
%     legend(rvNames)
%     
%     figure                      % convergence of SIGMA
%     for DimInp = 1:M
%         semilogx(N_seq,VarEEPlot(:,DimInp),'-x')
%         hold on
%     end
%     hold off
%     title('\sigma_{EE}')
%     legend(rvNames)
% end
% 
% %%
% % Print out on Command Window the total amount of time used for the
% % Sensitivity Assessment of the problem at issue
% timeT = toc(timeT);
% message = 'Total time: %4.2f seconds \n\n';
% fprintf(message,timeT)
% 
% %%
% % Elementary Effect Evaluation interfaced with COMSOL software
% %
% % Author: Gian Marco Melito         Date: December 2018
% % Project: Aorta.COMSOL
% % TU Graz