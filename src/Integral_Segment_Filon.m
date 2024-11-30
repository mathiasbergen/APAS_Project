classdef Integral_Segment_Filon
    % Class definition on generalized Filon integration
    
    properties (Access=private)
        f_Filon;
        g_Filon;
        Filon_Integrate = @(A,B,C,K,L,M,eta_sym)exp(M.*-1i-L.*eta_sym.*1i-K.*eta_sym.^2.*1i).*((B.*5.0e-1i)./K-A.*1.0./K.^2.*L.*2.5e-1i)+ ...
            (A.*eta_sym.*exp(M.*-1i-L.*eta_sym.*1i-K.*eta_sym.^2.*1i).*5.0e-1i)./K+1.0./K.^2.*sqrt(pi).*exp(M.*-1i+(L.^2.*2.5e-1i)./K).* ...
            erfi((L.*5.0e-1i+K.*eta_sym.*1i).*1.0./sqrt(K.*-1i)).*1.0./sqrt(K.*-1i).*(A.*L.^2.*1i+C.*K.^2.*4.0i+A.*K.*2.0-B.*K.*L.*2.0i).*1.25e-1i;
        n_coarse = 3;
        n_fine = 5;
    end
    properties (Access=public)
        C_INT;
    end
    
    methods
        function obj = Integral_Segment_Filon(f_Filon,g_Filon,C_INT)   % contructor
            obj.f_Filon=f_Filon;
            obj.g_Filon=g_Filon;
            obj.C_INT = C_INT;
        end
        
        function out = DoIntegrationCoarse(obj,eta,eta_step)
            eta_vec = linspace(eta,eta+eta_step,obj.n_coarse);
            
            FF = polyfit(eta_vec,obj.f_Filon(eta_vec),2);
            GG = polyfit(eta_vec,obj.g_Filon(eta_vec),2);
            out = obj.Filon_Integrate( FF(1,1), FF(1,2), FF(1,3), GG(1,1) , GG(1,2) , GG(1,3) , eta_vec(end)) ...
                -obj.Filon_Integrate( FF(1,1), FF(1,2), FF(1,3), GG(1,1) , GG(1,2) , GG(1,3) , eta_vec(1));
        end
        
        function out = DoIntegrationFine(obj,eta,eta_step)
            eta_vec = linspace(eta,eta+eta_step,obj.n_fine);
            n1 = ceil(obj.n_fine/2);
            eta_vec1 = eta_vec(1:n1);
            eta_vec2 = eta_vec(n1:obj.n_fine);
            
            FF1 = polyfit(eta_vec1,obj.f_Filon(eta_vec1),2);
            GG1 = polyfit(eta_vec1,obj.g_Filon(eta_vec1),2);
            FF2 = polyfit(eta_vec2,obj.f_Filon(eta_vec2),2);
            GG2 = polyfit(eta_vec2,obj.g_Filon(eta_vec2),2);
            %             (type: Ax^2+Bx+C)
            out1 = obj.Filon_Integrate( FF1(1,1), FF1(1,2), FF1(1,3), GG1(1,1) , GG1(1,2) , GG1(1,3) , eta_vec1(end)) ...
                -obj.Filon_Integrate( FF1(1,1), FF1(1,2), FF1(1,3), GG1(1,1) , GG1(1,2) , GG1(1,3) , eta_vec1(1));
            out2 = obj.Filon_Integrate( FF2(1,1), FF2(1,2), FF2(1,3), GG2(1,1) , GG2(1,2) , GG2(1,3) , eta_vec2(end)) ...
                -obj.Filon_Integrate( FF2(1,1), FF2(1,2), FF2(1,3), GG2(1,1) , GG2(1,2) , GG2(1,3) , eta_vec2(1));
            
            out = out1+out2;
        end
        
        function [IntegralSum,eta,etaStep] = adaptFilon(obj,etaStep,eta,etaEndSafe,errTresh)
            IntegralSum = 0;
            warning('')
% eta 
% etaEndSafe
% etaStep
            while eta<etaEndSafe
                if eta+etaStep>etaEndSafe
                    etaStep =  etaEndSafe - eta;
                end
%                 etaStep
%                 return
                Filon_segment_course = obj.DoIntegrationCoarse(eta, etaStep);
                Filon_segment_fine = obj.DoIntegrationFine(eta, etaStep);
%                 lastwarn
                %                     abs(Filon_segment_course - Filon_segment_fine)
                if ~isempty(lastwarn)
                    break
                end
                if isnan(Filon_segment_course) || isnan(Filon_segment_fine)
                    etaStep = etaStep*0.93745644864384684;
                    Filon_segment_course = obj.DoIntegrationCoarse(eta, etaStep);
                    Filon_segment_fine = obj.DoIntegrationFine(eta, etaStep);
                end
                
                if abs(Filon_segment_course - Filon_segment_fine)*obj.C_INT  <=  errTresh
                    IntegralSum = IntegralSum + Filon_segment_fine;
                    eta = eta+etaStep;
                    etaStep = etaStep*2;
                else
                    etaStep=etaStep/2;
                end
            end %eta
        end
        
    end
end



