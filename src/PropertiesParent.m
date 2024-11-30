classdef PropertiesParent
    
    properties (Constant)
        rho_f = 1000;   % Fluid density
        v0 = 1;         % Piston particle velocity
        K_f_real = 2.205225000000000e+09;  % Fluid real bulk modulus.  K_f_real  = real(c_f^2)  * rho_f
        mu_p_real = 7.837520000000000e+10; % Plate real shear modulus. mu_p_real = real(c_s^2) * rho_p
        K_p_real = 1.6276693333333333333333333333333333333e+11;  % Plate real bulk modulus.  K_p_real  = real(c_p^2) * rho_p - 4/3*mu_p_real
        rho_p = 8000;   % Plate density
        d = 6.05e-3;
        z1 = 270e-3;    % Plate distance from piston source
        loss_K_f = struct('model', 'Constant Q', 'Value',inf); % Parameter after 'Value' sets constant Q factor for the fluid bulk modulus.
        loss_K_p = struct('model', 'Constant Q', 'Value',inf); % Parameter after 'Value' sets constant Q factor for the plate bulk modulus.
        loss_mu_p = struct('model', 'Constant Q', 'Value',inf);% Parameter after 'Value' sets constant Q factor for the plate shear modulus.
    end
    
    properties (Access=public)
        k_p;
        h_p;
        h_f;
        K_f;
        mu_p;
        K_p;
        eta_complexFactor;
        f;
	end
    
    methods
        function obj = PropertiesParent(f)   % contructor
            if nargin == 1
                obj.f = f;
            else
                f = 0;
                obj.f = f;
            end
            obj.K_f = obj.K_f_real.*obj.set_complexFactor_K_f;
            obj.K_p = obj.K_p_real.*obj.set_complexFactor_K_p;
            obj.mu_p = obj.mu_p_real.*obj.set_complexFactor_mu_p;
            obj.eta_complexFactor = 1./sqrt(obj.set_complexFactor_K_f);
     
            c_f = sqrt(obj.K_f./obj.rho_f);
            c_ps = sqrt(obj.mu_p./obj.rho_p);
            c_pl = sqrt( (obj.K_p + 4/3*obj.mu_p)./obj.rho_p);
                        
            w = 2*pi*f;
            obj.k_p = w/c_ps;
            obj.h_p = w/c_pl;
            obj.h_f = w/c_f;
        end
        
        function complexFactor_K_f = set_complexFactor_K_f(obj)
            if strcmp(obj.loss_K_f.model,'Constant Q')
                Q = obj.loss_K_f.Value;
                complexFactor_K_f = 1+i/Q;
            elseif strcmp(obj.loss_K_f.model,'Varying Q')
                % implement method
            end
        end
        
        function complexFactor_K_p = set_complexFactor_K_p(obj)
            if strcmp(obj.loss_K_p.model,'Constant Q')
                Q = obj.loss_K_p.Value;
                complexFactor_K_p = 1+i/Q;
            elseif strcmp(obj.loss_K_p.model,'Varying Q')
                % implement method
            end
        end
        
        function complexFactor_mu_p = set_complexFactor_mu_p(obj)
            if strcmp(obj.loss_mu_p.model,'Constant Q')
                Q = obj.loss_mu_p.Value;
                complexFactor_mu_p = 1+i/Q;
            elseif strcmp(obj.loss_mu_p.model,'Varying Q')
                % implement method
            end
        end
    end
end
