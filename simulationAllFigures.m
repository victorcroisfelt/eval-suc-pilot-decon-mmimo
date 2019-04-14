%This Matlab script can be used to generate all figures in the article:
%
%Victor Croisfelt Rodrigues, Taufik Abrão, "An Evaluation of Successive 
%Pilot Decontamination in Massive MIMO", Semina: Ciências Exatas e 
%Tecnológicas, Londrina, v. 39, n. 2, ago./dez. 2018.
%
%Download paper: https://doi.org/10.5433/1679-0375.2018v39n2p107
%
%This is version 2.0 (Last edited: 04-14-2019)
%
%License: This code is licensed under the GPLv3 license. If you in any way
%use this code for research that results in publications, please reference 
%our original article as shown above.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Initialization
close all;
clearvars;

%% Simulation parameters

%Choose the desired simulation subfigure:
%   simulation == 1: to simulate Figure 1
%   simulation == 2: to simulate Figures 2 and 3
%
simulation = 1;

%Number of BSs
L = 7;

%Define parameters in each simulation scenario
if simulation == 1

    %Number of BS antennas
    M = 100;
    
    %Determine maximum number of BS antennas
    Mmax = M;

    %Select range of total uplink transmit power per UE [mW]
    rhoulRange = logspace(-3,3,10);
    
    %Number of UEs per BS
    K = 70;
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(rhoulRange);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 1000;
    
elseif simulation == 2
    
    %Select the range of number of BS antennas
    Mrange = round(logspace(2,5,10));
    
    %Determine maximum number of BS antennas
    Mmax = max(Mrange);
    
    %Define total uplink transmit power per UE [mW]
    rhoul = 200;
    
    %Number of UEs per BS
    K = 70;
    
    %Number of measurement points along horizontal axis in final figure
    nbrOfPoints = length(Mrange);
    
    %Number of setups with randomly generated statistics
    nbrOfSetups = 1000;
    
end

%% Scenario Setup

%Define BS radius [m]
cellRadius = 500;

%Distance between BSs [m]
interBSdistance = sqrt(3)*cellRadius;

%Define BS positions using complex coordinates [m]
BSlocations = [0+1i*0 interBSdistance*exp(1i*(pi/3*(0:5)))];
%Important: The center cell or home cell is considered to be index 1,
%i.e., j = 1.

%% Propagation parameters

%Pathloss exponent
alpha = 3.76;

%Standard deviation of large-scale fading variations (shadowing) [dB]
std_shad = 8;

%Communication bandwidth [Hz]
bandwidth = 20e6;

%Define noise figure at BS [dB]
noiseFigure = 10;

%Compute total noise power [dBm]
noiseVariancedBm = -174 + 10*log10(bandwidth) + noiseFigure;

%Select length of coherence block [samples]
tauc = 200;

%% Simulation

%Prepare to save the different simulation results
if simulation == 1
    
    %NMSE - Normalized Mean-Squared Error
    meanNMSE_class_sc = zeros(nbrOfPoints,nbrOfSetups); % classical single-cell
    meanNMSE_class = zeros(nbrOfPoints,nbrOfSetups); % classical 
    meanNMSE_eval1 = zeros(nbrOfPoints,nbrOfSetups); % evaluated 1
    meanNMSE_eval2 = zeros(nbrOfPoints,nbrOfSetups); % evaluated 2
   
elseif simulation == 2
    
    %SE - Spectral Efficiency
    meanSE_class_sc = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_class = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_eval1 = zeros(nbrOfPoints,nbrOfSetups);
    meanSE_eval2 = zeros(nbrOfPoints,nbrOfSetups);

    meanSumSE_class_sc = zeros(nbrOfPoints,nbrOfSetups);
    meanSumSE_class = zeros(nbrOfPoints,nbrOfSetups);
    meanSumSE_eval1 = zeros(nbrOfPoints,nbrOfSetups);
    meanSumSE_eval2 = zeros(nbrOfPoints,nbrOfSetups);
    
end

%Go through all setups
for s = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(s) ' setups out of ' num2str(nbrOfSetups)]);

    %Randomly distribute the UEs in the cell area
    UElocations = functionDistributeUniformlyUEs(L,K,cellRadius,BSlocations);
    
    %Prepare to store pathloss coefficients [dB]
    pathgaindB = zeros(K,L);
    
    %Go through all the cells
    for l = 1:L
       
        %Compute distances between UEs in BS l and BS j
        distancesBSj = abs(UElocations(:,l)-BSlocations(1));
        
        %Compute distant-dependent path gains [dB]
        pathgaindB(:,l) = - alpha*10*log10(distancesBSj);
        
    end
    
    %Compute the normalized channel gains, where the normalization is by
    %noise power
    channelGaindB = pathgaindB - noiseVariancedBm;
         
    %Compute the large-scale coefficients by adding the shadowing term
    betas = 10.^(std_shad.*randn(K,L)./10) .* 10.^(channelGaindB./10);
       
    %Go through all horizontal points
    for m = 1:nbrOfPoints

        %Output simulation progress
        disp([num2str(m) ' points out of ' num2str(nbrOfPoints)]);        
               
        if simulation == 1  
            
            %Extract current uplink transmitted power
            rhoul = rhoulRange(m);
           
        elseif simulation == 2
           
            %Extract current number of antennas
            M = Mrange(m);
                     
        end
        
        if simulation == 1
            
            %Compute the MMSE channel estimation
            [~,~,~,~,~,~,NMSE_class_sc,NMSE_class,NMSE_eval1,NMSE_eval2] = functionMMSEChannelEstimates(L,K,rhoul,betas);
            
            %Store the NMSE-based results
            meanNMSE_class_sc(m,s) = mean(NMSE_class_sc);
            meanNMSE_class(m,s) = mean(NMSE_class(:,1));
            meanNMSE_eval1(m,s) = mean(NMSE_eval1);
            meanNMSE_eval2(m,s) = mean(NMSE_eval2);

        elseif simulation == 2
            
            %Compute the MMSE channel estimation
            [taup,taupl,psis_class_sc,psis_class,psis_eval1,psis_eval2] = functionMMSEChannelEstimates(L,K,rhoul,betas);
            
            %Compute UL (uplink) SE
            [SE_class_sc,SE_class,SE_eval1,SE_eval2] = functionComputeSE_UL(L,M,K,tauc,taup,taupl,rhoul,betas,psis_class_sc,psis_class,psis_eval1,psis_eval2);
            
            %Store the UL SE per user
            meanSE_class_sc(m,s) = mean(SE_class_sc);
            meanSE_class(m,s) = mean(SE_class);
            meanSE_eval1(m,s) = mean(SE_eval1);
            meanSE_eval2(m,s) = mean(SE_eval2);
            
            %Store the UL sum SE per user
            meanSumSE_class_sc(m,s) = sum(SE_class_sc);
            meanSumSE_class(m,s) = sum(SE_class);
            meanSumSE_eval1(m,s) = sum(SE_eval1);
            meanSumSE_eval2(m,s) = sum(SE_eval2);
       
        end

    end
    
end

%% Plot simulation results

if simulation == 1
    
    %Plot Figure 1
    figure;
    hold on; box on;
    
    plot(rhoulRange,10*log10(mean(meanNMSE_class_sc,2)),'k:','LineWidth',1.5)
    plot(rhoulRange,10*log10(mean(meanNMSE_class,2)),'k--','LineWidth',1.5)
    plot(rhoulRange,10*log10(mean(meanNMSE_eval1,2)),'-.^','LineWidth',1)
    plot(rhoulRange,10*log10(mean(meanNMSE_eval2,2)),'-.v','LineWidth',1)
    
    ylabel('Average NMSE [dB/user/antenna]');
    xlabel('Radiated power per terminal [mW]');

    legend('Classical: single-cell','Classical (3)','Estimator 1 (13)','Estimator 2 (17)','Location','best');
    
    set(gca,'XScale','log');
    
elseif simulation == 2
    
    %Plot Figure 2
    figure;
    hold on; box on;
        
    plot(Mrange,mean(meanSE_class_sc,2),'k:','LineWidth',1.5)
    plot(Mrange,mean(meanSE_class,2),'k--','LineWidth',1.5)
    plot(Mrange,mean(meanSE_eval1,2),'-.^','LineWidth',1)
    plot(Mrange,mean(meanSE_eval2,2),'-.v','LineWidth',1)

    ylabel('Average UL SE [bit/s/Hz/user]');
    xlabel('Number of BS antennas ($M$)');
    
    legend('Classical: single-cell','Classical (3)','Estimator 1 (13)','Estimator 2 (17)','Location','best');
    
    set(gca,'XScale','log');
  
    %Plot Figure 3
    figure;
    subplot(1,1,1)
    hold on; box on;

    plot(Mrange,mean(meanSumSE_class_sc,2),'k:','LineWidth',1.5)
    plot(Mrange,mean(meanSumSE_class,2),'k--','LineWidth',1.5)
    plot(Mrange,mean(meanSumSE_eval1,2),'-.^','LineWidth',1)
    plot(Mrange,mean(meanSumSE_eval2,2),'-.v','LineWidth',1)
    
    ylabel('Average sum UL SE [bit/s/Hz]');
    xlabel('Number of BS antennas ($M$)');

    legend('Classical: single-cell','Classical (3)','Estimator 1 (13)','Estimator 2 (17)','Location','best');
    
    set(gca,'XScale','log');  
    
end