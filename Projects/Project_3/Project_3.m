% Gopal Gupta
% Roll Number - 20CH30008
clc;
clear;
warning off;

encoding_decoding;
function encoding_decoding
% loading the Data
load("data_cn_project_iii_a22.mat");
%% QUESTION 1: Auto Correlation
xr = xcorr(Stimulus,Stimulus,50);
figure(1);
plot(-50:50,xr/20000,'k');
hold on;
title('Autocorrelation function of the Stimulus');
xlabel('𝜏 (in ms)');
ylabel('R(𝜏)');
plot(-50:50,zeros(1,101));
set(gca,'FontSize',20);
hold off;

%% QUESTION 2: PSTH PLOTING
PSTH = zeros(4,20000);
train_psth = zeros(4,20);
for neuron_no = 1:4
    for trials = 1:50
        PSTH(neuron_no,:) = PSTH(neuron_no,:) + histcounts(All_Spike_Times{neuron_no,trials}*1000,0:20000)*1000;
    end
    PSTH(neuron_no,:) = PSTH(neuron_no,:)/50;
    for j = 1:20
        train_psth(neuron_no,j) =  train_psth(neuron_no,j) + mean(PSTH(neuron_no,1+(j-1)*1000:j*1000));
    end
end
% Plotting the PSTH
figure(2);
for neuron_no = 1:4
    subplot(4,1,neuron_no);
    plot(train_psth(neuron_no,1:20),'k');
    xlabel('Time (sec)');
    ylabel('R(t) (spikes/s)');
    title(['PSTHs or Mean firing rate of Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',13);
end
%% QUESTION 3: Poisson or non-Poisson distribution

bin_size = [10 20 50 100 200 500];
for bin = 1:length(bin_size)
    figure(bin+2);
    for neuron_no = 1:4
        mean_trial = zeros(4,50);
        variance_trial = zeros(4,50);
        for trials = 1:50
            spikes = histcounts(All_Spike_Times{neuron_no,trials}*1000,0:bin_size(bin):20000);
            mean_trial(neuron_no,trials) = mean(spikes);
            variance_trial(neuron_no,trials)=var(spikes);
        end
        ax__max_lm = max([max(mean_trial),max(variance_trial)]);
        ax__min_lm = min([min(mean_trial),min(variance_trial)]);
        subplot(2,2,neuron_no);
        scatter(variance_trial,mean_trial);
        hold on;
        fano_factor_line = [0,ax__max_lm];
        plot(fano_factor_line,fano_factor_line,'r');
        xlabel('Mean');ylabel('Variance');
        title(['Mean vs Variance plot Neuron ' num2str(neuron_no) ' with  bin size ' num2str(bin_size(bin)) ' ms' ]);
        xlim([ax__min_lm ax__max_lm]);
        ylim([ax__min_lm ax__max_lm]);
        set(gca,'FontSize',13);
    end
end

%% Question 4 : Spike Triggered Average
sta = zeros(4,100);
h_t = zeros(4,100);
% Calculating the correction factor
Rxx = autocorr(Stimulus,99);
%Rxx = xcorr(Stimulus,Stimulus,99)/20000;
Css = toeplitz(Rxx);
for neuron_no = 1:4
    spike_count = 0;
    for trials = 1:50
        spike_count = spike_count + nnz(All_Spike_Times{neuron_no,trials} <= 15);
        for spi_train = 1:nnz(All_Spike_Times{neuron_no,trials} <= 15)
            rnd_sp_time = round(All_Spike_Times{neuron_no,trials}(spi_train)*1000);        
            if rnd_sp_time < 100
                stm = Stimulus(1:rnd_sp_time);
                compl_stm = [zeros(1,100-length(stm)) stm];
                sta(neuron_no,:) = sta(neuron_no,:) + compl_stm;
            else
                sta(neuron_no,:) = sta(neuron_no,:) + Stimulus((rnd_sp_time-99):rnd_sp_time);
            end
        end
    end
    sta(neuron_no,:) = sta(neuron_no,:)/spike_count;
    figure(9);
    subplot(2,2,neuron_no);
    plot(0:99,fliplr(sta(neuron_no,:)),0:99,zeros(1,100));
    xlabel('Time (ms)');
    ylabel('STA');
    ylim([-0.2 0.2]);
    title(['Filter h(t) without correction factor of Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',13);
    % Applying Whitening Correction
    h_t(neuron_no,:) = fliplr((Css\sta(neuron_no,:)')');
    figure(10);
    subplot(2,2,neuron_no);
    plot(0:99,h_t(neuron_no,:),0:99,zeros(1,100));
    xlabel('Time (ms)');
    ylabel('STA');
    ylim([-1.3 1.3]);
    title(['Filter h(t) with correction factor of Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',13);
end    

%% Question 5: Determining the output non-linearity
% Evaluating the prediction using the filter
x_pred = zeros(4,15099);
for neuron_no = 1:4
    x_pred(neuron_no,:) = conv(Stimulus(1:15000),h_t(neuron_no,:));
end
% Taking the trained data for analysis
x_train = zeros(4,15000);
for neuron_no = 1:4
    x_train(neuron_no,:) = x_pred(neuron_no,1:15000);
end
% y_train Data set
y_train = zeros(4,15000);
for neuron_no = 1:4
y_train(neuron_no,:) = PSTH(neuron_no,1:15000);
end
% A mean plot of λ(t) for binned y(t) values
x_mean_test = zeros(4,300);
y_mean_test = zeros(4,300);
bin_x_test = 50;
bn_sz_test = 1:bin_x_test:15000;
for neuron_no = 1:4
    for j = 1:ceil(15000/bin_x_test) 
        x_mean_test(neuron_no,j) = mean(x_train(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
        y_mean_test(neuron_no,j) = mean(y_train(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
    end
    figure(11);
    subplot(2,2,neuron_no);
    scatter(x_mean_test(neuron_no,:),y_mean_test(neuron_no,:),5,'filled');
    xlabel('y(t)');
    ylabel('λ(t)');
    title(['Estimated y(t) vs Measured λ(t) for Neuron ' num2str(neuron_no)]);
    grid on
    set(gca,'FontSize',13);
end
% finding the nonlinear output function from y_mean and x_mean
[fit_const,fit_error_values] = fitting_data(x_mean_test,y_mean_test);
save('fits', 'fit_const', 'fit_error_values');


%% Question 6: Prediction performance and pruning of filter parameters

% Evaluating the prediction using the filter
% linear filter
x_test_pre = zeros(4,5000);
y_test = zeros(4,5000);
for neuron_no = 1:4
    pred = conv(Stimulus(15001:20000),h_t(neuron_no,:));
    x_test_pre(neuron_no,:) = pred(1:5000);
end
% Applying the non linear filter
for neuron_no = 1:4
    for i = 1:5000
        x_test_pre(neuron_no,i) = (fit_const{neuron_no,1}.a)/(1+exp(-fit_const{neuron_no,1}.b*(x_test_pre(neuron_no,i)-fit_const{neuron_no,1}.c)));
    end
end
% Getting the actual λ(t) 
for neuron_no = 1:4
        y_test(neuron_no,:) = PSTH(neuron_no,15001:20000);
end
%  Evaluating the r_squared value

r_square = zeros(4,1);
fprintf('The r_squared value before Pruning');
for neuron_no = 1:4
    r_value = corrcoef(y_test(neuron_no,:),x_test_pre(neuron_no,:));
    r_square(neuron_no) = r_value(2)^2;
    fprintf('The r_squared value of neurone %d is',neuron_no);
    disp(r_square(neuron_no));
    %figure(13);
    %subplot(2,2,neuron_no);
    %scatter(x_test(neuron_no,:),y_test(neuron_no,:),5,'blue','filled');
    %xlabel('Predicted λ(t)');
    %ylabel('Actual λ(t)');
    %title(['Predicted λ(t) vs Actual λ(t) for test set of Neuron ' num2str(neuron_no)]);
end
% Plotting the Estimated and Actual firing rate
x_mean_test_pre = zeros(4,100);
y_mean_test_pre = zeros(4,100);
bin_x_test = 50;
bn_sz_test = 1:bin_x_test:5000;
for neuron_no = 1:4
    for j = 1:ceil(5000/bin_x_test) 
        x_mean_test_pre(neuron_no,j) = mean(x_test_pre(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
        y_mean_test_pre(neuron_no,j) = mean(y_test(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
    end
    figure(13);
    subplot(2,2,neuron_no);
    plot(bn_sz_test,y_mean_test_pre(neuron_no,:));
    hold on;
    plot(bn_sz_test,x_mean_test_pre(neuron_no,:));
    ylabel('λ(t)');
    xlabel('Time (ms)');
    title(['Predicted λ(t) vs Actual λ(t) for test set of Neuron ' num2str(neuron_no)]);
    legend('Actual','Predicted');
    hold off;
    grid on;
    set(gca,'FontSize',13);
end
% Since r_square value of the 1 and 4 neuron is very low nad we discard  them
% pruning the h(2,:) and h(3,:)
counter = zeros(100,2);
cor_coef = zeros(100,2);
x_test_pruned = zeros(2,5000);
for neuron_no = 2:3
    initi_r_sqr = r_square(neuron_no);
    prev_r_sqr = r_square(neuron_no);
    next_r_sqr = prev_r_sqr;
    index = 0;
    count = 0;
    while (prev_r_sqr - next_r_sqr) < 0.014
        prev_r_sqr = next_r_sqr;
        min_va = max(h_t(neuron_no,:));
        for i = 1:100
            if abs(h_t(neuron_no,i))<min_va && abs(h_t(neuron_no,i)) ~= 0
                min_va = abs(h_t(neuron_no,i));
                index = i;
            end
        end
        if index == 0
            break;
        end
        h_t(neuron_no,index) = 0;
        index = 0;
        pred_pruned = conv(Stimulus(15001:20000),h_t(neuron_no,:));
        x_test_pruned(neuron_no-1,:) = pred_pruned(1:5000);
        for i = 1:5000
            x_test_pruned(neuron_no-1,i) = (fit_const{neuron_no,1}.a)/(1+exp(-fit_const{neuron_no,1}.b*(x_test_pruned(neuron_no-1,i)-fit_const{neuron_no,1}.c)));
        end
        r_value = corrcoef(y_test(neuron_no,:),x_test_pruned(neuron_no-1,:));
        next_r_sqr = r_value(2)^2;
        count = count + 1;
        counter(count,neuron_no-1) = count;
        cor_coef(count,neuron_no-1) = next_r_sqr;
    end
    figure(12+neuron_no);
    scatter(counter(:,neuron_no-1),cor_coef(:,neuron_no-1),'filled');
    ylabel('Prediction performance');
    xlabel('number of iteration');
    title(['Plot of Prediction performance vs no of iteration for Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',16);
 end

% Ploting the filter
four_trn = zeros(4,100);
four_axis = (-50:49)';
for neuron_no = 1:4  
    figure(16);
    subplot(2,2,neuron_no);
    plot(0:99,h_t(neuron_no,:),0:99,zeros(1,100));
    xlabel('Time (ms)');
    ylabel('h(t)');
    ylim([-1.3 1.3]);
    title(['Filter h(t) after pruning of Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',13);
end
% Plotting the FFT of the Filter
for neuron_no = 1:4
    y_shift = circshift(h_t(neuron_no,:),[0,length(h_t(neuron_no,:))/2]);
    four_trn(neuron_no,:) = fft(y_shift,length(h_t(neuron_no,:)))/length(h_t(neuron_no,:));
    Y_fft = circshift(four_trn(neuron_no,:),[0,length(h_t(neuron_no,:))/2]);
    figure(17);
    subplot(2,2,neuron_no);
    plot(four_axis,abs(Y_fft),four_axis,zeros(1,100));
    xlabel('freqeuncy');
    ylabel('h(f)');
    title(['Fourier transform of filter of Neuron ' num2str(neuron_no)]);
    set(gca,'FontSize',13);
end
%
x_test = zeros(4,5000);
for neuron_no = 1:4
    pred = conv(Stimulus(15001:20000),h_t(neuron_no,:));
    x_test(neuron_no,:) = pred(1:5000);
end
% Applying the non linear filter
for neuron_no = 1:4
    for i = 1:5000
        x_test(neuron_no,i) = (fit_const{neuron_no,1}.a)/(1+exp(-fit_const{neuron_no,1}.b*(x_test(neuron_no,i)-fit_const{neuron_no,1}.c)));
    end
end
fprintf('The r_squared value after Pruning');
for neuron_no = 1:4
    r_value = corrcoef(y_test(neuron_no,:),x_test(neuron_no,:));
    r_square(neuron_no) = r_value(2)^2;
    %disp(r_value);
    fprintf('The r_squared value of neurone %d is',neuron_no);
    disp(r_square(neuron_no));

    %figure(19);
    %subplot(2,2,neuron_no);
    %plot(y_test(neuron_no,:));
    %hold on;
    %plot(x_test(neuron_no,:));
    %xlabel('λ(t)');
    %ylabel('Time (ms)');
    %title(['Predicted λ(t) vs Actual λ(t) for test set of Neuron ' num2str(neuron_no)]);
    %legend('Actual','Predicted');
    %hold off;
end
%
x_mean_test = zeros(4,100);
y_mean_test = zeros(4,100);
bin_x_test = 50;
bn_sz_test = 1:bin_x_test:5000;
for neuron_no = 1:4
    for j = 1:ceil(5000/bin_x_test) 
        x_mean_test(neuron_no,j) = mean(x_test(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
        y_mean_test(neuron_no,j) = mean(y_test(neuron_no,bn_sz_test(j):bn_sz_test(j)+(bin_x_test-1)));
    end
    figure(18);
    subplot(2,2,neuron_no);
    plot(bn_sz_test,y_mean_test(neuron_no,:));
    hold on;
    plot(bn_sz_test,x_mean_test(neuron_no,:));
    ylabel('λ(t)');
    xlabel('Time (ms)');
    title(['Predicted λ(t) vs Actual λ(t) for test set of Neuron ' num2str(neuron_no)]);
    legend('Actual','Predicted');
    hold off;
    grid on;
    set(gca,'FontSize',13);
end

%% Question 7B : Discrimination based on Victor and Purpura (VP) Spike Distance Metric (SDM) 
% Different value of the cost 
q = [0, 0.001, 0.01, 0.1, 1, 10, 100];
mut_inf = zeros(4,100,length(q));
% Running the operaion for 100 times
for sp_case = 1:10
    % Randomly Choosing 8 response/spike_train_segment
    rand_time = randperm(19901,8);
    for neuron_no = 1:4
        spike_time_segment = cell(8,50);
        for trials = 1:50
            spike_train = All_Spike_Times{neuron_no,trials};
            for rnd_spike = 1:8
                spike_time_segment{rnd_spike,trials} = spike_train(spike_train >= rand_time(rnd_spike)/1000 & spike_train < (rand_time(rnd_spike)+100)/1000);
            end
        end
        
        for cost_val = 1:length(q)
            conf_matr = zeros(8, 8);
            for rnd_spike = 1:8
                for trials = 1:50
                    mean_SD = mean_spike_distance(rnd_spike,trials,spike_time_segment,q(cost_val));
                    mean_SD(rnd_spike) = mean_SD(rnd_spike)*50/49;
                    [~,ind] = min(mean_SD);
                    conf_matr(rnd_spike,ind) = conf_matr(rnd_spike,ind) + 1;
                end
            end
            conf_matr = conf_matr/50;
            mut_inf(neuron_no,sp_case,cost_val) = mut_inf(neuron_no,sp_case,cost_val) + Mut_Info(conf_matr)/10;
        end    
    end
end
%Calculating the Confidence Interval
CI90 = CI(0.9,mut_inf);
a = mut_inf(:,1,:);
b = CI90(:,1,:);
% Plotting the Descrimination 
figure(19);
for neuron_no = 1:4
    subplot(2,2,neuron_no);
    plot(log10(q), a(neuron_no,:));
    hold on;
    plot(log10(q), a(neuron_no,:)-b(neuron_no,:), 'r--');
    plot(log10(q), a(neuron_no,:)+b(neuron_no,:), 'r--');
    [~,p] = max(a(neuron_no,:));
    p = plot(log10(q(p)), a(neuron_no,1,p), '+');
    set(p, 'linewidth', 2);
    hold off;
    title(['Discrimination - Neuron ' num2str(neuron_no)]);
    xlabel('log_{10}(q)');
    ylabel('Mean MI(q) (90% confidence intervals)');
    legend('Mean MI','Confidence Interval','Location','northwest');
end
end
% Finding the mean VP SDM
function mean_SP = mean_spike_distance(rnd_spike,trials,spike_time_segment,q)
    mean_SP = zeros(1,8);
    for rnd_spike_1 = 1:8
        for trials_1 = 1:50
            if (trials_1 == trials && rnd_spike_1 == rnd_spike)
                continue
            end
            mean_SP(rnd_spike_1) = mean_SP(rnd_spike_1) + spkd(spike_time_segment{rnd_spike,trials},spike_time_segment{rnd_spike_1,trials_1},q);
        end
    end
    mean_SP = mean_SP/50;
end
% Function to find the Mutual Information of the Confusion matrix
function MI = Mut_Info(conf_matr)
    MI = 0;
    for i = 1:size(conf_matr,1)
        for j = 1:size(conf_matr,2)
            if (conf_matr(i,j) ~= 0)
                MI = MI + conf_matr(i,j)/size(conf_matr,1)*log2(conf_matr(i,j)/sum(conf_matr(:,j))); % confusion matrix has entries of p(y/x)
            end
        end
    end
end
% Function to find the confidence Interval
function CI90 = CI(interval,mut_inf)
    %MIbar = mean(mut_inf,2);
    MIstd = std(abs(mut_inf),0,2);
    alpha = 1 - interval;
    T_multiplier = tinv(1-alpha/2, 99);
    CI90 = T_multiplier*MIstd/sqrt(99);
end
% Function To find the Spike Distace Metric
function d=spkd(tli,tlj,cost)
%
% d=spkd(tli,tlj,cost) calculates the "spike time" distance
% (Victor & Purpura 1996) for a single cost
%
% tli: vector of spike times for first spike train
% tlj: vector of spike times for second spike train
% cost: cost per unit time to move a spike
%
%  Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
%  Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
%
nspi=length(tli);
nspj=length(tlj);

if cost==0
   d=abs(nspi-nspj);
   return
elseif cost==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
%
%     INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
%
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if (nspi && nspj)
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+cost*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end

% finding the nonlinear output function from y_mean and x_mean fit function
function [fit_const,fit_error_values] = fitting_data(x_mean,y_mean)
    fit_const = cell(4,1);
    fit_error_values = struct( 'sse', cell( 4, 1 ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    % Fitting for Neuron 1
    [x_data,y_data] = prepareCurveData(x_mean(1,:),y_mean(1,:));
    % Using the sigmoid curve to fit and setting the fittype options
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.7513 0.2551 0.5060];
    % Fit model to data.
    [fit_const{1}, fit_error_values(1)] = fit( x_data, y_data, ft, opts);
    % Plot fit with data.
    figure(12);
    subplot(2,2,1);
    plot( fit_const{1}, x_data, y_data );
    title('Fit for neuron 1');
    xlabel('y(t)');
    ylabel('λ(t)');
    grid on;
    set(gca,'FontSize',13);

    % Fitting for Neuron 2
    [x_data,y_data] = prepareCurveData(x_mean(2,:),y_mean(2,:));
    % Using the sigmoid curve to fit and setting the fittype options
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.3500 0.1966 0.2511];
    % Fit model to data.
    [fit_const{2}, fit_error_values(2)] = fit( x_data, y_data, ft, opts);
    % Plot fit with data.
    subplot(2,2,2);
    plot( fit_const{2}, x_data, y_data);
    title('Fit for neuron 2');
    xlabel('y(t)');
    ylabel('λ(t)');
    grid on;
    set(gca,'FontSize',13);

    % Fitting for Neuron 3
    [x_data,y_data] = prepareCurveData(x_mean(3,:),y_mean(3,:));
    % Using the sigmoid curve to fit and setting the fittype options
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.7537 0.3804 0.3804];
    % Fit model to data.
    [fit_const{3}, fit_error_values(3)] = fit( x_data, y_data, ft, opts);
    % Plot fit with data.
    subplot(2,2,3);
    plot( fit_const{3}, x_data, y_data );
    title('Fit for neuron 3');
    xlabel('y(t)');
    ylabel('λ(t)');
    grid on;
    set(gca,'FontSize',13);

    % Fitting for Neuron 4
    [x_data,y_data] = prepareCurveData(x_mean(4,:),y_mean(4,:));
    % Using the sigmoid curve to fit and setting the fittype options
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.3112 0.5285 0.1656];
    % Fit model to data.
    [fit_const{4}, fit_error_values(4)] = fit( x_data, y_data, ft, opts);
    % Plot fit with data.
    subplot(2,2,4);
    plot( fit_const{4}, x_data, y_data );
    title('Fit for neuron 4');
    xlabel('y(t)');
    ylabel('λ(t)');
    grid on;
    set(gca,'FontSize',13);
end





