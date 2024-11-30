% APAS main program.
% Sets a few input parameters and loads more prameters in PropertiesParent.m

% Programmed by Mathias Myrtveit Sæther
% (C) 2024 Mathias Myrtveit Sæther. This file is free software;  you can redistribute   
% it and/or modify it only under the the terms of the MIT LICENSE which must be included along with this file. 


%% Program starts
clear
warning('') % Clear lastwarn. A MATLAB lastwarn warning on badly scaled polynomial fitting is expected when using Filon method.
% The warning is OK because it triggers MATLAB to execute Gauss method instead of Filon.
warning('off','all')

%% Input parameters
f_vec  = 600e3;         % Frequency
z_vec = 376.05e-3;           % Distance from receiver to source
r_vec = 0;             % Lateral distance; or radius of receiver when receiverType = 'transducer'
a_vec = 10.55e-3;         % Radius of source
recAndModelType = 'hydrophone';          % Input argument: 'transducer'  or  'hydrophone' or 'diffCorrTransducer' , 'diffCorrHydrophone' , 'planeWave' , 'KinslerAndFrey'
wavePropSetup = 'transm';   % Input argument: 'transm','transm1','transm2','echo','echo1','echo2','freefield',
%% Input parameters numerical integration scheme
errTresh = 1e-1;       % Absolute local error tolerance in the Gauss and the Filon-type adaptive methods except around the integrand peak. Typically set to 0.01
etaStepBreak = 1e-8;  % Search for singuarities when eta step is less than etaStepBreak. Typical value 1e-8
etaStepInit = 1e-2;   % Initial staring step size
integrationBoundaryDeltaAbove_h_f = 200;
filonMethodStarts_OscPeriodCutOff_or_EtaCutoff = 0.1; % apparant period. Typical value: 0.1-0.2. Defines when Filon method kicks in. High value, Filon Method starts early
% For a value set in the range <0,10>, a eta period in the integrand is defined when Filon method kicks in. For a value set
% outside <0,10>, Filon method kicks in at eta for the set value (for inf, the Gauss is always used, for 0, Filon is always used)

%% Loop over field
if (strcmp(recAndModelType,'transducer') || strcmp(recAndModelType,'diffCorrTransducer')) && r_vec == 0
    disp('Error:: Transducer radius is not defined!!!')
    return
end

c_wa0 = sqrt(PropertiesParent.K_f_real/PropertiesParent.rho_f);
plateVarsSave = struct(PropertiesParent);

scale = 1;
PlaneWaveOrKF_model = 0;
etaStep  = etaStepInit;
eta = 0;
tic

for hh = 1:length(a_vec)
    a = a_vec(hh);
    for ii = 1:length(r_vec)
        r = r_vec(ii);
        for jj = 1:length(z_vec)
            z = z_vec(jj);
            for kk = 1:length(f_vec)
                f = f_vec(kk);
                f
                eta = 0;
                if strcmp(recAndModelType,'transducer') || strcmp(recAndModelType,'diffCorrTransducer')
                    integrandFunctions = TransducerFunctions(a,r,z,f);
                    if strcmp(wavePropSetup,'transm')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm;
                        f_Filon = @integrandFunctions.f_FilonTransm;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        planeWave = PlaneWavePropagation(f);
                        TR_coeff_0 = integrandFunctions.TransmissionCoeff(0);
                    elseif strcmp(wavePropSetup,'transm1')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm1;
                        f_Filon = @integrandFunctions.f_FilonTransm1;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        TR_coeff_0 = integrandFunctions.TransmissionCoeffNoRefl(0);
                    elseif strcmp(wavePropSetup,'transm2')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm2;
                        f_Filon = @integrandFunctions.f_FilonTransm2;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        TR_coeff_0 = integrandFunctions.TransmissionCoeffNoRefl2(0);
                    elseif strcmp(wavePropSetup,'echo')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho;
                        integrand_term = @integrandFunctions.termPartsGaussEcho;
                        f_Filon = @integrandFunctions.f_FilonEcho;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeff(0);
                    elseif strcmp(wavePropSetup,'echo1')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho1;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho1;
                        integrand_term = @integrandFunctions.termPartsGaussEcho1;
                        f_Filon = @integrandFunctions.f_FilonEcho1;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeffNoRefl(0);
                    elseif strcmp(wavePropSetup,'echo2')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho2;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho2;
                        integrand_term = @integrandFunctions.termPartsGaussEcho2;
                        f_Filon = @integrandFunctions.f_FilonEcho2;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeffNoRefl(0);
                    elseif strcmp(wavePropSetup,'freefield')
                        integrand_org = @integrandFunctions.integrandGaussFreefield;
                        integrand_parts = @integrandFunctions.integrandPartsGaussFreefield;
                        integrand_term = @integrandFunctions.termPartsGaussFreefield;
                        f_Filon = @integrandFunctions.f_FilonFreefield;
                        g_Filon = @integrandFunctions.g_FilonFreeField;
                        z_fluid = z;
                        TR_coeff_0 = 1;
                    end
                    
                    if strcmp(recAndModelType,'diffCorrTransducer')
                        planeWave = PlaneWavePropagation(f);
                        scale = planeWave.planeWaveFluid(z_fluid).* TR_coeff_0;
                    end

                elseif strcmp(recAndModelType,'hydrophone') || strcmp(recAndModelType,'diffCorrHydrophone')
                    integrandFunctions = HydrophoneFunctions(a,r,z,f);
                    if strcmp(wavePropSetup,'transm')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm;
                        f_Filon = @integrandFunctions.f_FilonTransm;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        TR_coeff_0 = integrandFunctions.TransmissionCoeff(0);                        
                    elseif strcmp(wavePropSetup,'transm1')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm1;
                        f_Filon = @integrandFunctions.f_FilonTransm1;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        TR_coeff_0 = integrandFunctions.TransmissionCoeffNoRefl(0);                        
                    elseif strcmp(wavePropSetup,'transm2')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussTransm2;
                        f_Filon = @integrandFunctions.f_FilonTransm2;
                        g_Filon = @integrandFunctions.g_FilonTransmAll;
                        z_fluid = z-PropertiesParent.d;
                        TR_coeff_0 = integrandFunctions.TransmissionCoeffNoRefl2(0);                        
                    elseif strcmp(wavePropSetup,'echo')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho;
                        integrand_term = @integrandFunctions.termPartsGaussEcho;
                        f_Filon = @integrandFunctions.f_FilonEcho;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeff(0);                        
                    elseif strcmp(wavePropSetup,'echo1')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho1;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho1;
                        integrand_term = @integrandFunctions.termPartsGaussEcho1;
                        f_Filon = @integrandFunctions.f_FilonEcho1;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeffNoRefl(0);                        
                    elseif strcmp(wavePropSetup,'echo2')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        integrand_org = @integrandFunctions.integrandGaussEcho2;
                        integrand_parts = @integrandFunctions.integrandPartsGaussEcho1;
                        integrand_term = @integrandFunctions.termPartsGaussEcho1;
                        f_Filon = @integrandFunctions.f_FilonEcho2;
                        g_Filon = @integrandFunctions.g_FilonEchoAll;
                        z_fluid = PropertiesParent.z1*2-z;
                        TR_coeff_0 = integrandFunctions.ReflectionCoeffNoRefl(0);                        
                    elseif strcmp(wavePropSetup,'freefield')
                        integrand_org = @integrandFunctions.integrandGaussFreefield;
                        integrand_parts = @integrandFunctions.integrandPartsGaussFreefield;
                        integrand_term = @integrandFunctions.termPartsGaussFreefield;
                        f_Filon = @integrandFunctions.f_FilonFreefield;
                        g_Filon = @integrandFunctions.g_FilonFreeField;
                        z_fluid = z;
                        TR_coeff_0 = 1;                        
                    end
                    
                    if strcmp(recAndModelType,'diffCorrHydrophone')
                        planeWave = PlaneWavePropagation(f);
                        scale = planeWave.planeWaveFluid(z_fluid).* TR_coeff_0;
                    end
                    
                    
                elseif strcmp(recAndModelType,'planeWave')
                    integrand_org = [];
                    f_Filon = [];
                    g_Filon = [];
                    eta = inf;
                    integrandFunctions = HydrophoneFunctions(a,r,z,f);
                    planeWave = PlaneWavePropagation(f);
                    if strcmp(wavePropSetup,'transm')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        z_fluid = z-PropertiesParent.d;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.TransmissionCoeff(0);
                    elseif strcmp(wavePropSetup,'transm1')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        z_fluid = z-PropertiesParent.d;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.TransmissionCoeffNoRefl(0);
                    elseif strcmp(wavePropSetup,'transm2')
                        if z < PropertiesParent.z1; disp('Error!!! z is smaller than z1!!!'); return; end
                        z_fluid = z-PropertiesParent.d;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.TransmissionCoeffNoRefl2(0);
                    elseif strcmp(wavePropSetup,'echo')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        z_fluid = PropertiesParent.z1*2-z;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.ReflectionCoeff(0);
                    elseif strcmp(wavePropSetup,'echo1')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        z_fluid = PropertiesParent.z1*2-z;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.ReflectionCoeffNoRefl(0);
                    elseif strcmp(wavePropSetup,'echo2')
                        if z > PropertiesParent.z1; disp('Error!!! z is larger than z1!!!'); return; end
                        z_fluid = PropertiesParent.z1*2-z;
                        PW = planeWave.planeWaveFluid(z_fluid).*planeWave.ReflectionCoeffNoRefl2(0);
                    elseif strcmp(wavePropSetup,'freefield')
                        z_fluid = z;
                        PW = planeWave.planeWaveFluid(z_fluid);
                    end
                    
                    PlaneWaveOrKF_model = PW;
                    scale = a.*PropertiesParent.rho_f.*f.*2.*pi.*PropertiesParent.v0;
                    
                elseif strcmp(recAndModelType,'KinslerAndFrey')
                    integrand_org = [];
                    f_Filon = [];
                    g_Filon = [];
                    eta = inf;
                    z_fluid = z;
                    integrandFunctions = TransducerFunctions(a,r,z,f);
                    h_f = integrandFunctions.h_f;
                    
                    KF = PropertiesParent.rho_f *2*pi*f/h_f*  PropertiesParent.v0*  ( exp(-i*h_f.*z)-exp(-i*h_f.*sqrt(z.^2+a.^2) ) );
                    PlaneWaveOrKF_model = KF;
                    scale = a.*PropertiesParent.rho_f.*f.*2.*pi.*PropertiesParent.v0;
                end
                
                %% settings
                h_f0 = real(integrandFunctions.h_f);
                h_p0 = real(integrandFunctions.h_p);
                k_p0 = real(integrandFunctions.k_p);
                etaEnd = h_f0+integrationBoundaryDeltaAbove_h_f;
                
                eta_cutOf_Shampine = sqrt(  h_f0^2  -  ( 2*pi/z_fluid + sqrt(h_f0^2-[0:1:h_f0].^2)  ).^2 ) - [0:1:h_f0]; % find appearant period in fluid
                indx_tmp = find(imag(eta_cutOf_Shampine) == 0,1);  % find when imaginary part is zero
                eta_cutoff = find(abs(eta_cutOf_Shampine(indx_tmp:end))<filonMethodStarts_OscPeriodCutOff_or_EtaCutoff,1)+indx_tmp-1;
                etaStartFilon = eta_cutoff;
                if isempty(eta_cutoff)
                    eta_cutoff = h_f0;
                end
                if filonMethodStarts_OscPeriodCutOff_or_EtaCutoff >= 10 || filonMethodStarts_OscPeriodCutOff_or_EtaCutoff == 0
                    eta_cutoff = filonMethodStarts_OscPeriodCutOff_or_EtaCutoff;
                end
                if eta_cutoff > h_f0    % needed if eta_cutoff is manually set
                    eta_cutoff = h_f0;
                end
                %% settings
                C_INT = a.*PropertiesParent.rho_f.*f.*2.*pi.*PropertiesParent.v0;
                integralGauss = Integral_Segment_Gauss(integrand_org,C_INT);
                integralFilon = Integral_Segment_Filon(f_Filon,g_Filon,C_INT);
                if exist('integrand_parts','var')
                    integralParts = Integral_Segment_Gauss(integrand_parts,C_INT);
                end
                
                %% Adaptive Gauss integration
                IntegralSum1=0;
                IntegralSumFilon=0;
                IntegralSum2 = 0;
                IntegralSumSafeEvanescent=0;
                IntegralSumEvanescent = [];
                
                [IntegralSum1,eta,etaStep]  = integralGauss.adaptIntegration(etaStepInit,eta,eta_cutoff,errTresh,etaStepBreak);
                if eta_cutoff ~= h_f0
                    [IntegralSumFilon,eta,etaStep] = integralFilon.adaptFilon(etaStep,eta,h_f0-etaStepBreak,errTresh);
                end
                [IntegralSum2,eta,etaStep]  = integralGauss.adaptIntegration(etaStep,eta,h_f0,errTresh,etaStepBreak);
                
                %% safe val evanescent
                if etaStep<etaStepBreak & ~strcmp(wavePropSetup,'transm') & ~strcmp(wavePropSetup,'transm1') & ~strcmp(wavePropSetup,'transm2')
                    delta_eta = h_f0-eta;
                    eta1 = eta;
                    eta2 = eta1+3*delta_eta;
                    [IntegralSum3_Gauss,eta,etaStep]  = integralParts.adaptIntegration(etaStep,eta1,eta2,errTresh,0);
                    IntegralSum3_term = integrandFunctions.termPartsGaussEcho(eta2)-integrandFunctions.termPartsGaussEcho(eta1);
                    IntegralSumSafeEvanescent = IntegralSum3_term-IntegralSum3_Gauss;
                end
                safeValEvanescent = eta;
                safeValEvanescent0 = safeValEvanescent;
                etaStep0_Evanescent = etaStepBreak;  % TJA
                
                %%  Org
                while safeValEvanescent<etaEnd
                    while true
                        %% Find 1st peak in evanescent region
                        [IntegralSum_Ev_All_tmp,eta_peakCrash,etaStep]  = integralGauss.adaptIntegration(etaStep0_Evanescent,safeValEvanescent,etaEnd,errTresh,etaStepBreak);
                        if eta_peakCrash == etaEnd
                            IntegralSumEvanescent = [IntegralSumEvanescent,IntegralSum_Ev_All_tmp];
                            safeValEvanescent = etaEnd;
                            break
                        else
                            infoStruct = peakDetectSub(eta_peakCrash,etaStepBreak,integrand_org,safeValEvanescent,etaEnd);
                            if infoStruct.isPeak == 0
                                safeValEvanescent=etaEnd;
                                break
                            end
                            etaPeak1 = infoStruct.etaPeak;
                            delta_eta_peak1 = etaPeak1-safeValEvanescent;
                            eta1_peak1 = safeValEvanescent;
                            eta2_peak1 = etaPeak1+delta_eta_peak1;
                            integrand_dual_peak1 = @(eta) integrand_org(eta) + integrand_org(safeValEvanescent+safeValEvanescent-eta+2*delta_eta_peak1);
                            integralGauss_dual_peak1 = Integral_Segment_Gauss(integrand_dual_peak1,C_INT);
                        end
                        
                        %% Find 2nd peak in evanescent region
                        [IntegralSum_Ev_peak1_tmp,eta_peakCrash,etaStep]  = integralGauss_dual_peak1.adaptIntegration(etaStep0_Evanescent,safeValEvanescent,etaPeak1,errTresh,etaStepBreak);
                        if eta_peakCrash == etaPeak1         % Far out
                            IntegralSumEvanescent = [IntegralSumEvanescent,IntegralSum_Ev_peak1_tmp];
                            safeValEvanescent = eta2_peak1;
                            break
                        else                                 % flipped
                            infoStruct = peakDetectSub(eta_peakCrash,etaStepBreak,integrand_dual_peak1,safeValEvanescent,etaEnd);
                            etaPeak2 = infoStruct.etaPeak;
                            delta_eta_peak2_tmp1 = etaPeak1-etaPeak2;
                            delta_eta_peak2_tmp2 = etaPeak2-safeValEvanescent;
                            delta_eta_peak2 = min([delta_eta_peak2_tmp1,delta_eta_peak2_tmp2]);
                            eta1_peak2 = etaPeak2-delta_eta_peak2;
                            eta2_peak2 = etaPeak2+delta_eta_peak2;
                            integrand_dual_peak2 = @(eta) integrand_dual_peak1(eta) + integrand_dual_peak1(eta1_peak2+eta1_peak2-eta+2*delta_eta_peak2);
                            integralGauss_dual_peak2 = Integral_Segment_Gauss(integrand_dual_peak2,C_INT);
                        end
                        
                        [IntegralSum_Ev_peak2a,eta_peakCrash,etaStep]  = integralGauss_dual_peak2.adaptIntegration(etaStep0_Evanescent,eta1_peak2,etaPeak2,errTresh,etaStepBreak);
                        if delta_eta_peak2_tmp1 < delta_eta_peak2_tmp2
                            eta1_rest = safeValEvanescent0;
                            eta2_rest = eta1_peak2;
                        else
                            eta1_rest = eta2_peak2;
                            eta2_rest = etaPeak1;
                        end
                        [IntegralSum_Ev_peak2b,eta_peakCrash,etaStep]  = integralGauss_dual_peak1.adaptIntegration(etaStep0_Evanescent,eta1_rest,eta2_rest,errTresh,etaStepBreak);
                        
                        IntegralSumEvanescent = [IntegralSumEvanescent,IntegralSum_Ev_peak2a,IntegralSum_Ev_peak2b];
                        safeValEvanescent = eta2_peak1;
                        break
                        
                    end
                end
                
                %%
                eta = safeValEvanescent;
                
                %% End evanescent
                IntegralSum_total = PlaneWaveOrKF_model + IntegralSum1 + IntegralSumFilon + IntegralSum2 + IntegralSumSafeEvanescent + sum(IntegralSumEvanescent);
                tot_press = IntegralSum_total*C_INT;
                freq_spectrum{hh,ii,jj,kk} = tot_press/scale;
                
            end % f
        end % z
    end % r
end % a
compute_time = toc;

% save APAS_result.mat freq_spectrum a_vec f_vec z_vec r_vec errTresh compute_time recAndModelType wavePropSetup plateVarsSave
resultsFolder = fullfile(pwd, '..', 'results');
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end
resultFilePath = fullfile(resultsFolder, 'APAS_result.mat');
save(resultFilePath, 'a_vec', 'f_vec', 'z_vec', 'r_vec', 'errTresh', 'compute_time', 'recAndModelType', 'wavePropSetup', 'plateVarsSave');

%% Quick access of first simulation (comment out)
pressure = squeeze(cell2mat(freq_spectrum(1,1,1,1)));
20*log10(abs(pressure))
% abs(pressure)
% angle(pressure)


