function [diff_data] = timeCal_modelDataDifference(time_vec,output_matrix,fittedData,titleTXT,N)
% compute and plot the difference between the desired model output and the
% fitted data coming from the MRI data


% Normalise model time
time_vec_norm = (time_vec - min(time_vec))./(max(time_vec) - min(time_vec));
time_array = [1 11 21 31 41 51 61]; % indices for model instances

diff_data.inTime = zeros(length(time_vec_norm),size(output_matrix,2));
% Difference MRI - Model
for idx = 1:size(output_matrix,2)
    diff_data.inTime(:,idx) = (output_matrix(:,idx) - fittedData)./abs(fittedData).*100;
    diff_data.inTime(1,idx) = 0;
    diff_data.difference(:,idx) = (output_matrix(:,idx) - fittedData);
end

% Get distribution of relative error in time
diff_data.q25_data_inTime = quantile(diff_data.inTime,0.25,2); % q25
diff_data.q50_data_inTime = quantile(diff_data.inTime,0.50,2); % q50 = median
diff_data.q75_data_inTime = quantile(diff_data.inTime,0.75,2); % q75

% Find one number to evaluate difference
% diff_data.indicator = mean(diff_data.inTime);
% diff_data.indicatorSSE = sum(diff_data.difference(time_array,:).^2);
% likelihood function: higher exponent N -> higher weights to good runs;
% L = 1/( (1/(2*obs))*sum_t(model(t) - data(t))^2 )^N
% diff_data.likelihood = sum(diff_data.difference(time_array,:))./(2*7).^N;
% diff_data.likelihood = (var(diff_data.difference(time_array,:),[],1)).^(-N);
% rescale likelihood
% diff_data.likelihood = diff_data.likelihood./sum(diff_data.likelihood);

% Test formula agreement
diff_data.likelihood = (1./((sum(diff_data.difference,1)).^2)./(2*61)).^N;
% rescale likelihood
diff_data.likelihood = diff_data.likelihood./sum(diff_data.likelihood);


%% Plot
figure('Visible','off')

% MRI and mean model response
subplot(3,2,[1,2])
plot(linspace(0,1,length(fittedData)),fittedData,'k','LineWidth',1)
hold on
plot(time_vec_norm,mean(output_matrix,2),'b','LineWidth',1)
plot(time_vec_norm,median(output_matrix,2),'r','LineWidth',1)
ylabel(titleTXT)
legend('MRI','model mean','median','Location','best')
xlabel('$t/T$ [-]')
grid on
if strcmpi(titleTXT(1:2),'SA')
    set(gca, 'YScale', 'log');
elseif strcmpi(titleTXT(1:3),'Vol')
    set(gca, 'YScale', 'log');
end

% relative error in time
subplot(3,2,[3,4])
plot(time_vec_norm,diff_data.q50_data_inTime,'r','LineWidth',1)
hold on
plot(time_vec_norm,mean(diff_data.inTime,2),'b','LineWidth',1)
plot(time_vec_norm,diff_data.q25_data_inTime,':r','LineWidth',1)
plot(time_vec_norm,diff_data.q75_data_inTime,':r','LineWidth',1)
ylabel('$\varepsilon$ [\%]')
legend('Median','Mean','$q_{25}$; $q_{75}$','Location','best')
xlabel('$t^*_\mathrm{mod}$ [-]')
grid on

% pdf of likelihood
subplot(325)
histogram(log(diff_data.likelihood),'Normalization','probability')
% set(gca,'XScale','log')
ylabel('$f(L)$')
xlabel('$log(L)$')
set(gca,'YScale','log')
grid on

% box-plot of indicator
subplot(326)
boxplot(diff_data.likelihood)
ylabel('$log(L)$')
set(gca,'YScale','log')
grid on



end