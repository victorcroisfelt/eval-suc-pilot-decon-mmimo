function [taup,taupl,psis_class_sc,psis_class,psis_eval1,psis_eval2,NMSE_class_sc,NMSE_class,NMSE_eval1,NMSE_eval2] = functionMMSEChannelEstimates(L,K,rhoul,betas)
%Compute the normalized mean-squared error (NMSE) for the classical and
%multi pilot training phases approaches based on the closed form
%expressions given by the minimum mean square error (MMSE) channel
%estimator. 'psis' correspods to the channel estimate average variance.
%
%This Matlab function is used in the article:
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
%@Inputs:
%   L: Number of BSs or cells.
%   K: Number of UEs per cell.
%   rhoul: Uplink transmit power per UE (same for every user).
%   betas: K x L matrix with the average large-scale coefficient of the
%   users over the entire system in relation to the center cell.
%
%@Outputs:
%   taup: Length of the pilots.
%   taupl: Length of the pilots for the multi pilot training scenario.
%   psis_class_sc: K x 1 vector with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell
%   for the classical scheme considering only cell j.
%   psis_class: K x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell 
%   for the classical scheme.
%   psis_eval1: K/L x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell
%   for the evaluated scheme 1 (multiple pilot training phases).
%   psis_eval2: K/L x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell
%   for the evaluated scheme 2 (multiple pilot training phases).
%   NMSE_class_sc: K x 1 vector with the NMSE of each user in relation to
%   cell j for the classical scheme.
%   NMSE_class: K x L matrix with the NMSE of each user in relation to cell
%   j for the classical scheme considering only cell j.
%   NMSE_eval1: K/L x L matrix with the NMSE of each user in relation to
%   cell j for the evaluated scheme 1 (multiple pilot training phases).
%   NMSE_eval2: K/L x L matrix with the NMSE of each user in relation to
%   cell j for the evaluated scheme 2 (multiple pilot training phases).
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%% Classical schemes

%Length of pilot sequences [samples]
taup = K;

%Compute the normalized pilot transmitted power [mW]
rho_p = taup*rhoul;

%Prepare to store estimated variances
psis_class_sc = zeros(K,1);
psis_class = zeros(K,L);

if nargout > 6
    
    %Prepare to store the NMSE for each user within each cell
    NMSE_class_sc = zeros(K,1);
    NMSE_class = zeros(K,L);
    
end

%Go through all UEs in cell j
for k = 1:K
    
    %Compute the denominator
    denom_class = 1+rho_p*sum(betas(k,:));
    denom_class_sc = 1+rho_p*betas(k,1);
    
    %Compute the estimated variance
    psis_class_sc(k) = (rho_p*betas(k,1)^2)/denom_class_sc;
    
    if nargout > 6
        
        %Compute the NMSE per user per antenna
        NMSE_class_sc(k) = (betas(k,1)-psis_class_sc(k))/betas(k,1);
        
    end
    
    %Go through all cells
    for l = 1:L
        
        %Compute the estimated variance
        psis_class(k,l) = (rho_p*betas(k,l)^2)/denom_class;
        
        if nargout > 6
            
            %Compute the NMSE per user per antenna
            NMSE_class(k,l) = (betas(k,l)-psis_class(k,l))/betas(k,l);
            
        end
        
    end
    
end

%% Evaluated schemes

%Length of pilot sequences [samples]
taupl = K/L;

%Compute the normalized pilot transmitted power [mW]
rho_p = taupl*rhoul;

%Prepare to store estimation covariance matrix
psis_eval1 = zeros(K/L,1);
psis_eval2 = zeros(K/L,1);

if nargout > 6
    
    %Prepare to store the NMSE for each user within the center cell
    NMSE_eval1 = zeros(K/L,1);
    NMSE_eval2 = zeros(K/L,1);
    
end

%Go through all UEs considering multiple pilot training phases
for k = 1:K/L
    
    %Compute the denominator
    denom_eval1 = 2+rho_p*betas(k,1);
    denom_eval2 = (L^2-L+1)+L^2*rho_p*betas(k,1);
    
    %Compute the estimate variance
    psis_eval1(k) = (rho_p*betas(k,1)^2)/denom_eval1;
    psis_eval2(k) = (L^2*rho_p*betas(k,1)^2)/denom_eval2;
    
    if nargout > 6
        
        %Compute the NMSE per user per antenna
        NMSE_eval1(k) = (betas(k,1)-psis_eval1(k))/betas(k,1);
        NMSE_eval2(k) = (betas(k,1)-psis_eval2(k))/betas(k,1);
        
    end
    
end
