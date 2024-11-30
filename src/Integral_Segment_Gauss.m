classdef Integral_Segment_Gauss
    % Class definition on Gauss integration
    
    
    properties (Access=private)
        x_Gauss_points4 = [-0.861136311594053  -0.339981043584856   0.339981043584856   0.861136311594053];
        w_Gauss_weights4 = [0.347854845137454   0.652145154862546   0.652145154862546   0.347854845137454];
        x_Gauss_points5 = [-0.906179845938664  -0.538469310105683   0.000000000000000   0.538469310105683   0.906179845938664];
        w_Gauss_weights5 = [0.236926885056189   0.478628670499367   0.568888888888889   0.478628670499367   0.236926885056189];
        x_Gauss_points6 = [-0.9324695142031520278123, -0.661209386466264513661,-0.2386191860831969086305,0.238619186083196908631,0.661209386466264513661,0.9324695142031520278123];
        w_Gauss_weights6 = [0.1713244923791703450403,0.3607615730481386075698,0.4679139345726910473899,0.46791393457269104739,0.3607615730481386075698,0.1713244923791703450403];
        integrand;
    end
    properties (Access=public)
        C_INT;
    end
    
    methods
        function obj = Integral_Segment_Gauss(integrand,C_INT)   % contructor
            obj.integrand = integrand;
            obj.C_INT = C_INT;
        end
        
        function segment = Gauss4(obj,eta,eta_step)
            a = eta;
            b = eta+eta_step;
            eta_t = (a+b)./2+(b-a)./2.*obj.x_Gauss_points4;
            segment = (b-a)./2 .* obj.integrand(eta_t) * transpose(obj.w_Gauss_weights4);
        end
        
        function segment = Gauss5(obj,eta,eta_step)
            a = eta;
            b = eta+eta_step;
            eta_t = (a+b)./2+(b-a)./2.*obj.x_Gauss_points5;
            segment = (b-a)./2 .* obj.integrand(eta_t) * transpose(obj.w_Gauss_weights5);
        end
        
        function segment = Gauss6(obj,eta,eta_step)
            a = eta;
            b = eta+eta_step;
            eta_t = (a+b)./2+(b-a)./2.*obj.x_Gauss_points6;
            segment = (b-a)./2 .* obj.integrand(eta_t) * transpose(obj.w_Gauss_weights6);
        end
        
        function [IntegralSum,eta,etaStep] = adaptIntegration(obj,etaStep,eta,etaEnd,errTresh,etaBreak)
            IntegralSum = 0;
            while eta<etaEnd
                if eta+etaStep>etaEnd
                    etaStep =  etaEnd - eta;
                end
                if etaStep < etaBreak  % && historic?
                    break
                end
                gaussOrder4 = obj.Gauss4(eta,etaStep);
                gaussOrder5 = obj.Gauss5(eta,etaStep);
                
                %                 gaussOrder4
                %                 gaussOrder5
                %                                 eta
                %                                 etaStep
                %                 tt = 1
                
                if isnan(gaussOrder4) || isnan(gaussOrder5)
                    etaStep = etaStep*0.93745644864384684;
                    gaussOrder4 = obj.Gauss4(eta,etaStep);
                    gaussOrder5 = obj.Gauss5(eta,etaStep);
                end
                
                if abs( (gaussOrder4 - gaussOrder5)*obj.C_INT)  <=  errTresh
                    IntegralSum = IntegralSum + gaussOrder5;
                    eta = eta+etaStep;
                    etaStep = etaStep*2;
                    %                     tt1 = imag(obj.integrand( eta+etaStep*(0:0.2:1))  );
                    %                     tt2 = imag(obj.integrand( eta+etaStep*(0.2:0.2:1.2))  );
                    %                     (tt2-tt1)
                    
                else
                    etaStep=etaStep/2;
                end
            end
        end
        
        
        
        function [IntegralSum,eta0,eta_peakEnd,etaStep_peak] = adaptIntegrationPeak3(obj,etaBreak,errTresh,etaPeak)
            deltaEta_oneSide = etaBreak;
            eta = etaPeak - deltaEta_oneSide;
            etaStep_peak = 2*deltaEta_oneSide;
            gaussOrder4_0 = obj.Gauss4(eta,etaStep_peak);
            gaussOrder6_0 = obj.Gauss6(eta,etaStep_peak);
            
            while true
                %                 eta = eta - etaStep_peak/2;
                %             etaStep_peak = etaStep_peak*2;
                
                deltaEta_oneSide = deltaEta_oneSide*2;
                eta = etaPeak - deltaEta_oneSide;
                etaStep_peak = deltaEta_oneSide  *2;
                
                gaussOrder4 = obj.Gauss4(eta,etaStep_peak);
                gaussOrder6 = obj.Gauss6(eta,etaStep_peak);
                cond_sign = (sign(imag(gaussOrder4_0)) * sign(imag(gaussOrder4))*sign(imag(gaussOrder6_0))*sign(imag(gaussOrder6))) == 1;
                
                %                 obj.C_INT
%                 abs( (gaussOrder4 - gaussOrder6)*obj.C_INT)
%                 cond_sign
%                 deltaEta_oneSide
                
                
                if cond_sign && abs( (gaussOrder4 - gaussOrder6)*obj.C_INT)  <=  errTresh  && (abs(gaussOrder4) > abs(gaussOrder4_0) ) && (abs(gaussOrder6) > abs(gaussOrder6_0))
                    IntegralSum = gaussOrder6;
                    break
                else
                    gaussOrder4_0 = gaussOrder4;
                    gaussOrder6_0 = gaussOrder6;
                end
            end
            eta0 = eta;
            eta_peakEnd = eta+etaStep_peak;
            
        end
        
        
        
        function [IntegralSum,eta0,eta_peakEnd,etaStep_peak] = adaptIntegrationPeak3_parts(obj,eta,etaBreak,errTresh,etaPeak,integrand_term)
            etaStep0 = etaBreak; 
            etaStep = etaBreak;
            safeVal = eta;
            eta1 = etaPeak-etaStep*2;
            eta2 = etaPeak + etaStep;
            
%             gaussOrder4_0 = obj.Gauss4(eta1,etaStep) + obj.Gauss4(eta2,etaStep);
%                 gaussOrder6_0 = obj.Gauss6(eta1,etaStep) + obj.Gauss6(eta2,etaStep);
                gaussOrder4_0 = obj.Gauss4(eta1,etaBreak) + obj.Gauss4(eta2,etaBreak);
                gaussOrder6_0 = obj.Gauss6(eta1,etaBreak) + obj.Gauss6(eta2,etaBreak);
%                 gaussOrder4_0 = obj.Gauss4(eta1,etaStep0) + obj.Gauss4(eta2,etaStep0);
%                 gaussOrder6_0 = obj.Gauss6(eta1,etaStep0) + obj.Gauss6(eta2,etaStep0);
                
            while true
                etaStep0 = etaStep0*1.1;
            etaStep = etaStep*2;
            eta1 = eta1-etaStep;
%             eta1 = eta1-1e-9
            eta2 = eta2 + etaStep/2;

            %%

            
            %%
%                 gaussOrder4 = obj.Gauss4(eta1,etaStep) + obj.Gauss4(eta2,etaStep);
%                 gaussOrder6 = obj.Gauss6(eta1,etaStep) + obj.Gauss6(eta2,etaStep);
                gaussOrder4 = obj.Gauss4(eta1,etaBreak) + obj.Gauss4(eta2,etaBreak);
                gaussOrder6 = obj.Gauss6(eta1,etaBreak) + obj.Gauss6(eta2,etaBreak);
%                 gaussOrder4 = obj.Gauss4(eta1,etaStep0) + obj.Gauss4(eta2,etaStep0);
%                 gaussOrder6 = obj.Gauss6(eta1,etaStep0) + obj.Gauss6(eta2,etaStep0);
                
                cond_sign = (sign(imag(gaussOrder4_0)) * sign(imag(gaussOrder4))*sign(imag(gaussOrder6_0))*sign(imag(gaussOrder6))) == 1;
%                 abs( (gaussOrder4 - gaussOrder6)*obj.C_INT)
%                 eta1-safeVal
%                 gaussOrder4_parts = integrand_term(eta1+etaStep)-integrand_term(eta1) + integrand_term(eta2+etaStep)-integrand_term(eta2) - gaussOrder4
                gaussOrder4_parts = integrand_term(eta1+etaBreak)-integrand_term(eta1) + integrand_term(eta2+etaBreak)-integrand_term(eta2) - gaussOrder4;
%                 gaussOrder4_parts*obj.C_INT
                
                
%                 if abs( (gaussOrder4 - gaussOrder6)*obj.C_INT)  <=  errTresh
                if cond_sign && abs( (gaussOrder4 - gaussOrder6)*obj.C_INT)  <=  errTresh  && (abs(gaussOrder4) > abs(gaussOrder4_0) ) && (abs(gaussOrder6) > abs(gaussOrder6_0))
                    IntegralSum = gaussOrder6;
                    
                    [IntegralSum1,eta_out1,etaStep_out1] = adaptIntegration(obj,etaBreak,safeVal,eta1,errTresh,etaBreak);
                    term1 = integrand_term(eta1)-integrand_term(safeVal);
                    IntegralSum1_final = term1-IntegralSum1;
                    
                    [IntegralSum2,eta_out2,etaStep_out2] = adaptIntegration(obj,etaBreak,eta2+etaBreak,400,errTresh,etaBreak);
                    term2 = integrand_term(400)-integrand_term(eta2+etaBreak);
                    IntegralSum1_final2 = term2-IntegralSum2;
                    
                    break
                else
                    gaussOrder4_0 = gaussOrder4;
                    gaussOrder6_0 = gaussOrder6;
                end
            end
            eta0 = eta;
            eta_peakEnd = eta+etaStep_peak;
            
        end
        
    end
end
