function infoStruct = peakDetectSub(eta,etaStep,integrand,eta0,etaEnd)
%% Info
% Accurately founding a peak of type 1/x

%% Test
% % integrand = integrand_org;
% integrand = integrand_parts;
% % % eta = etaPeak;
% % % eta = 2116.1706;
% % eta = 1.481314559973524e+03;
% etaStep = 1e-8;
% % eta = 1.936632534319868e+03;
% etaEnd = h_f0+1000;
% eta0 = h_f0;

%% Prog starts
% eta0 = eta;
% eta = eta+etaStep;
if isnan(integrand(eta))
    eta = eta+2*etaStep;
end
integrand_vec = [];
etaStep_tmp = etaStep;
etaStepNano = etaStep;
isPeak = false;
% signFunInit = sign(imag(integrand(eta0)));
% signDerInit = sign( imag(integrand(eta0+etaStepNano)-integrand(eta0)));
% eta >= eta0  &&  eta<=etaEnd
while eta >= eta0  &&  eta<=etaEnd
    etaVec = [eta-etaStep,eta,eta+etaStep];
    integrand_vec = imag(integrand(etaVec));
%     etaStep
%     etaVec
%     integrand_vec
    if (integrand_vec(1)==integrand_vec(2)) && (integrand_vec(1) == integrand_vec(3))
        if etaStep == 0
            etaStep = etaStep_tmp;
        end
%         etaStep = etaStep*3;
        etaStep = etaStepNano;
        integrand_neg2 = (imag(integrand(eta-2*etaStep)));
        integrand_neg = (imag(integrand(eta-etaStep)));
        integrand_mid =  (imag(integrand(eta)));
        integrand_pos = (imag(integrand(eta+etaStep )));
        integrand_pos2 =(imag(integrand(eta+2*etaStep )));
        integrand_vec = [integrand_neg2 integrand_neg integrand_mid integrand_pos integrand_pos2];
        mid = ceil(length(integrand_vec)/2);
%         integrand_vec
        % peak 1/x type condition
        maxval = max(integrand_vec(:));
        idx_max_camdidates = find(integrand_vec == maxval);
        [a b] = min( abs(idx_max_camdidates-mid));
        mid_max = idx_max_camdidates(b);
        
        minval = min(integrand_vec(:));
        idx_min_camdidates = find(integrand_vec == minval);
        [a b] = min( abs(idx_min_camdidates-mid));
        mid_min = idx_min_camdidates(b);
        
        cond1 = abs(mid_max-mid_min) == 1;   %         1. side by side
        cond2 = sign(integrand_vec(mid_max)) ~= sign(integrand_vec(mid_min));  %         2. oposite sign
        cond3 = ~sum(ismember([mid_max mid_min],[1 5]));  %         3. in middle
        cond4 = ( sign(integrand_vec(1))==sign(integrand_vec(2)) && sign(integrand_vec(4))==sign(integrand_vec(5)) ); %         4. two outer elements have same sign
        
        if cond1 && cond2 && cond3 && cond4
            isPeak = true;
        end
        
        break
        
    elseif max(abs(integrand_vec)) == abs(integrand_vec(1))
        eta = eta-etaStep;
        etaStep = etaStep*2;
%         etaStep = etaStep*1.5;
    elseif max(abs(integrand_vec)) == abs(integrand_vec(3))
        eta = eta+etaStep;
        etaStep = etaStep*2;
%         etaStep = etaStep*1.5;
%         gg = 1
    elseif max(abs(integrand_vec)) == abs(integrand_vec(2))
        etaStep_tmp = etaStep;
        etaStep = etaStep/3;
    end
    
end


% eta_sub2 = eta
infoStruct.etaPeak = eta;
infoStruct.isPeak = isPeak;
infoStruct.etaStepPeakDetect = etaStep;
infoStruct.integrand_vec = integrand_vec;

