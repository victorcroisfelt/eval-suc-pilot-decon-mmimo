function [SE_class_sc,SE_class,SE_eval1,SE_eval2] = functionComputeSE_UL(L,M,K,tauc,taup,taupl,rhoul,betas,psis_class_sc,psis_class,psis_eval1,psis_eval2)
%Computation of the uplink spectral efficiency (SE) considering the
%evaluated methods presented throughout the paper.
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
%   M: Number of BS antennas.
%   K: Number of UEs per cell.
%   tauc: Length of coherence block.
%   taup: Length of the pilots.
%   taupl: Length of the pilots for the multi pilot training scenario.
%   rhoul: Uplink transmit power per UE (same for everyone).
%   betas: K x L matrix with the average large-scale coefficient of the 
%   users over the entire system in relation to the center cell. 
%   psis_class_sc: K x 1 vector with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell 
%   for the classical scheme considering only cell j. 
%   psis_class: K x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell
%   for the classical scheme.  
%   psis_eval1: K/L x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell 
%   for the proposed scheme 1 (multiple pilot training phases).
%   psis_eval2: K/L x L matrix with the average variance of the users'
%   channel estimate over the entire system in relation to the center cell
%   for the proposed scheme 2 (multiple pilot training phases).
%
%@Outputs:
%   SE_class_sc: K x 1 vector where element (k,1) is the uplink SE of UE k 
%   in cell j achieved with MR combining for the classical estimation 
%   approach considering only the operation of BS j.
%   SE_class: K x 1 vector where element (k,1) is the uplink SE of UE k 
%   in cell j achieved with MR combining for the classical estimation 
%   approach.
%   SE_eval1: K/L x 1 vector where element (k,1) is the uplink SE of UE  
%   k in cell j achieved with MR combining for the #1 proposed estimation 
%   approach considering only the operation of BS j.
%   SE_eval2: K/L x 1 vector where element (k,1) is the uplink SE of UE  
%   k in cell j achieved with MR combining for the #2 proposed estimation
%   approach considering only the operation of BS j.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

    %% Classical schemes

    %Compute the pre-log factor 
    prelogFactor_class = (tauc-taup)/tauc;

    %Prepare to store simulation results
    SE_class_sc = zeros(K,1);
    SE_class = zeros(K,1);

    %Go through all UEs
    for k = 1:K
        
        %Classical: single-cell
        signal = M*rhoul*psis_class_sc(k);
        interfnoise = 1+rhoul*sum(betas(k,1));
        
        SE_class_sc(k) = prelogFactor_class*real(log2(1+signal/interfnoise));

        %Classical: multi-cell
        signal = M*rhoul*psis_class(k,1);
        interfnoise = 1+rhoul*sum(betas(k,:))+M*rhoul*sum(psis_class(k,2:L));
      
        SE_class(k) = prelogFactor_class*real(log2(1+signal/interfnoise));
       
    end
    
    %% Evaluated schemes 

    %Compute the pre-log factor 
    prelogFactor_eval = (tauc-taupl)/tauc;

    %Prepare to store simulation results
    SE_eval1 = zeros(K/L,1);
    SE_eval2 = zeros(K/L,1);
    
    %Go through all UEs
    for k = 1:K/L
        
        %Evaluated 1
        signal = M*rhoul*psis_eval1(k,1);
        interfnoise = 1+rhoul*sum(betas(k,:))+M*rhoul*sum(psis_class(k,2:L));
      
        SE_eval1(k) = prelogFactor_eval*real(log2(1+signal/interfnoise));
       
     	%Evaluated 2  
        signal = M*rhoul*psis_eval2(k,1);
        interfnoise = 1+rhoul*sum(betas(k,:))+M*rhoul*sum(psis_class(k,2:L));
      
        SE_eval2(k) = prelogFactor_eval*real(log2(1+signal/interfnoise));
        
    end
   
end
