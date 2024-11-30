classdef TransmAndReflCoeff < PropertiesParent
    
    properties (Access=public)
        % From PropertiesParent
    end
    
    methods
        function obj = TransmAndReflCoeff(f)   % contructor
            obj@PropertiesParent(f); % invoking the superclass constructor
        end
        function T = TransmissionCoeff(obj,eta)
            h_pz = (-i*(-( obj.h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( obj.k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( obj.h_f.^2-eta.^2)).^(1/2));  % Compact form
            Y = obj.rho_f.*h_pz./(obj.rho_p.*h_fz).*obj.k_p.^4;
            S = (2.*eta.^2-obj.k_p.^2).^2.*cot(h_pz*obj.d/2)+4*eta.^2.*h_pz.*k_pz.*cot(k_pz.*obj.d/2);
            A = (2.*eta.^2-obj.k_p.^2).^2.*tan(h_pz*obj.d/2)+4*eta.^2.*h_pz.*k_pz.*tan(k_pz.*obj.d/2);
            T = -i*Y.*(A+S)./((S+i*Y).*(A-i*Y));      % convention: e^{+iwt}
        end
        
        function T = TransmissionCoeffNoRefl(obj,eta)
            rho_p = obj.rho_p;
            rho_f = obj.rho_f;
            k_p = obj.k_p;
            h_p = obj.h_p;
            h_f = obj.h_f;
            h_pz = (-i*(-( obj.h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( obj.k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( obj.h_f.^2-eta.^2)).^(1/2));  % Compact form
            d = obj.d;
            T = h_fz.*h_p.^2.*h_pz.*k_p.^2.*rho_f.*rho_p.*exp(d.*h_pz.*-1i).*(eta.^2-k_pz.^2).*(eta.^2+k_pz.^2).*(eta.^2.*h_p.^2.*-2.0+eta.^2.*k_p.^2+h_pz.^2.*k_p.^2).*1.0./(eta.^4.*h_fz.*h_p.^2.*rho_p.*2.0-eta.^4.*h_fz.*k_p.^2.*rho_p-eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*2.0-eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f+eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p+h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f+eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*4.0).^2.*-4.0;
        end
        
        function T = TransmissionCoeffNoRefl2(obj,eta)
            rho_p = obj.rho_p;
            rho_f = obj.rho_f;
            k_p = obj.k_p;
            h_p = obj.h_p;
            h_f = obj.h_f;
            h_pz = (-i*(-( obj.h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( obj.k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( obj.h_f.^2-eta.^2)).^(1/2));  % Compact form
            d = obj.d;
            T = h_fz.*h_p.^2.*h_pz.*k_p.^2.*rho_f.*rho_p.*exp(-d.*(h_pz.*3.0i+k_pz.*1i)).*(eta.^4-k_pz.^4).*(eta.^2.*h_p.^2.*-2.0+eta.^2.*k_p.^2+h_pz.^2.*k_p.^2).*1.0./(eta.^4.*h_fz.*h_p.^2.*rho_p.*2.0-eta.^4.*h_fz.*k_p.^2.*rho_p-eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*2.0-eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f+eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p+h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f+eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*4.0).^4.*(eta.^8.*h_fz.^2.*h_p.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^8.*h_fz.^2.*k_p.^4.*rho_p.^2.*exp(d.*k_pz.*1i)-eta.^8.*h_fz.^2.*h_p.^2.*k_p.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^4.*h_fz.^2.*h_p.^4.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^4.*h_fz.^2.*h_pz.^4.*k_p.^4.*rho_p.^2.*exp(d.*k_pz.*1i)-eta.^6.*h_fz.^2.*h_p.^4.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^6.*h_fz.^2.*h_pz.^2.*k_p.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*2.0+eta.^4.*h_p.^4.*h_pz.^2.*k_p.^4.*rho_f.^2.*exp(d.*k_pz.*1i)+eta.^4.*h_fz.^2.*k_p.^4.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i)-eta.^6.*h_fz.^2.*k_p.^4.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*2.0+h_fz.^2.*h_pz.^4.*k_p.^4.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i)+h_p.^4.*h_pz.^2.*k_p.^4.*k_pz.^4.*rho_f.^2.*exp(d.*k_pz.*1i)-eta.^4.*h_fz.^2.*h_p.^4.*h_pz.*k_pz.^3.*rho_p.^2.*exp(d.*h_pz.*1i).*3.2e+1+eta.^4.*h_fz.^2.*h_p.^4.*h_pz.*k_pz.^3.*rho_p.^2.*exp(d.*k_pz.*1i).*1.6e+1-eta.^6.*h_fz.^2.*h_p.^2.*h_pz.^2.*k_p.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^4.*h_fz.^2.*h_p.^4.*h_pz.^2.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*1.6e+1-eta.^4.*h_fz.^2.*h_p.^2.*k_p.^2.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^6.*h_fz.^2.*h_p.^2.*k_p.^2.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^2.*h_fz.^2.*h_pz.^2.*k_p.^4.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*2.0-eta.^2.*h_fz.^2.*h_pz.^4.*k_p.^4.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*2.0-eta.^4.*h_fz.^2.*h_pz.^2.*k_p.^4.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0+eta.^2.*h_p.^4.*h_pz.^2.*k_p.^4.*k_pz.^2.*rho_f.^2.*exp(d.*k_pz.*1i).*2.0+eta.^6.*h_fz.^2.*h_p.^4.*h_pz.*k_pz.*rho_p.^2.*exp(d.*h_pz.*1i).*3.2e+1-eta.^6.*h_fz.^2.*h_p.^4.*h_pz.*k_pz.*rho_p.^2.*exp(d.*k_pz.*1i).*1.6e+1+eta.^6.*h_fz.*h_p.^2.*h_pz.*k_p.^4.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*2.0-eta.^6.*h_fz.*h_p.^4.*h_pz.*k_p.^2.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*4.0+eta.^4.*h_fz.^2.*h_p.^2.*h_pz.*k_p.^2.*k_pz.^3.*rho_p.^2.*exp(d.*h_pz.*1i).*1.6e+1-eta.^4.*h_fz.^2.*h_p.^2.*h_pz.^3.*k_p.^2.*k_pz.*rho_p.^2.*exp(d.*h_pz.*1i).*1.6e+1-eta.^4.*h_fz.^2.*h_p.^2.*h_pz.*k_p.^2.*k_pz.^3.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^4.*h_fz.^2.*h_p.^2.*h_pz.^3.*k_p.^2.*k_pz.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^4.*h_fz.*h_p.^2.*h_pz.^3.*k_p.^4.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*2.0-h_fz.*h_p.^2.*h_pz.^3.*k_p.^4.*k_pz.^4.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*2.0+eta.^2.*h_fz.^2.*h_p.^2.*h_pz.^3.*k_p.^2.*k_pz.^3.*rho_p.^2.*exp(d.*h_pz.*1i).*1.6e+1-eta.^2.*h_fz.^2.*h_p.^2.*h_pz.^2.*k_p.^2.*k_pz.^4.*rho_p.^2.*exp(d.*k_pz.*1i).*4.0-eta.^2.*h_fz.^2.*h_p.^2.*h_pz.^3.*k_p.^2.*k_pz.^3.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^4.*h_fz.^2.*h_p.^2.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0-eta.^6.*h_fz.^2.*h_p.^2.*h_pz.*k_p.^2.*k_pz.*rho_p.^2.*exp(d.*h_pz.*1i).*1.6e+1+eta.^6.*h_fz.^2.*h_p.^2.*h_pz.*k_p.^2.*k_pz.*rho_p.^2.*exp(d.*k_pz.*1i).*8.0+eta.^2.*h_fz.*h_p.^4.*h_pz.^2.*k_p.^2.*k_pz.^3.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*8.0-eta.^2.*h_fz.*h_p.^2.*h_pz.*k_p.^4.*k_pz.^4.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*2.0+eta.^2.*h_fz.*h_p.^4.*h_pz.*k_p.^2.*k_pz.^4.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*4.0+eta.^4.*h_fz.*h_p.^4.*h_pz.^2.*k_p.^2.*k_pz.*rho_f.*rho_p.*exp(d.*k_pz.*1i).*8.0).*-4.0;
        end
        
        function R = ReflectionCoeff(obj,eta)
            h_pz = (-i*(-( obj.h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( obj.k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( obj.h_f.^2-eta.^2)).^(1/2));  % Compact form
            Y = obj.rho_f.*h_pz./(obj.rho_p.*h_fz).*obj.k_p.^4;
            S = (2.*eta.^2-obj.k_p.^2).^2.*cot(h_pz*obj.d/2)+4*eta.^2.*h_pz.*k_pz.*cot(k_pz.*obj.d/2);
            A = (2.*eta.^2-obj.k_p.^2).^2.*tan(h_pz*obj.d/2)+4*eta.^2.*h_pz.*k_pz.*tan(k_pz.*obj.d/2);
            R = (A.*S-Y.^2)./((S+i*Y).*(A-i*Y));      % convention: e^{+iwt}
        end
        
        function R = ReflectionCoeffNoRefl(obj,eta)
            rho_p = obj.rho_p;
            rho_f = obj.rho_f;
            k_p = obj.k_p;
            h_p = obj.h_p;
            h_f = obj.h_f;
            d = obj.d;
            h_pz = (-i*(-( h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( h_f.^2-eta.^2)).^(1/2));  % Compact form
            R = -(eta.^4.*h_fz.*h_p.^2.*rho_p.*-2.0+eta.^4.*h_fz.*k_p.^2.*rho_p+eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*2.0+eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f-eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p-h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f-eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*4.0)./(eta.^4.*h_fz.*h_p.^2.*rho_p.*2.0-eta.^4.*h_fz.*k_p.^2.*rho_p-eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*2.0-eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f+eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p+h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f+eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*4.0);
        end
        
        function R = ReflectionCoeffNoRefl2(obj,eta)
            rho_p = obj.rho_p;
            rho_f = obj.rho_f;
            k_p = obj.k_p;
            h_p = obj.h_p;
            h_f = obj.h_f;
            d = obj.d;
            h_pz = (-i*(-( h_p.^2-eta.^2)).^(1/2));  % Compact form
            k_pz = (-i*(-( k_p.^2-eta.^2)).^(1/2));  % Compact form
            h_fz = (-i*(-( h_f.^2-eta.^2)).^(1/2));  % Compact form
            R = h_fz.*h_p.^2.*h_pz.*k_p.^2.*rho_f.*rho_p.*exp(-d.*(h_pz.*2.0i+k_pz.*1i)).*(eta.^4-k_pz.^4).*(eta.^2.*h_p.^2.*-2.0+eta.^2.*k_p.^2+h_pz.^2.*k_p.^2).*1.0./(eta.^4.*h_fz.*h_p.^2.*rho_p.*2.0-eta.^4.*h_fz.*k_p.^2.*rho_p-eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*2.0-eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f+eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p+h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f+eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*4.0).^3.*(eta.^4.*h_fz.*h_p.^2.*rho_p.*exp(d.*k_pz.*1i).*-2.0+eta.^4.*h_fz.*k_p.^2.*rho_p.*exp(d.*k_pz.*1i)+eta.^2.*h_fz.*h_p.^2.*k_pz.^2.*rho_p.*exp(d.*k_pz.*1i).*2.0+eta.^2.*h_fz.*h_pz.^2.*k_p.^2.*rho_p.*exp(d.*k_pz.*1i)+eta.^2.*h_p.^2.*h_pz.*k_p.^2.*rho_f.*exp(d.*k_pz.*1i)-eta.^2.*h_fz.*k_p.^2.*k_pz.^2.*rho_p.*exp(d.*k_pz.*1i)-h_fz.*h_pz.^2.*k_p.^2.*k_pz.^2.*rho_p.*exp(d.*k_pz.*1i)+h_p.^2.*h_pz.*k_p.^2.*k_pz.^2.*rho_f.*exp(d.*k_pz.*1i)-eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*exp(d.*h_pz.*1i).*8.0+eta.^2.*h_fz.*h_p.^2.*h_pz.*k_pz.*rho_p.*exp(d.*k_pz.*1i).*4.0).*-4.0;
        end
    end
end



