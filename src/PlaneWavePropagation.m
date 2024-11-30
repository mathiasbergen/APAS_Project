classdef PlaneWavePropagation < TransmAndReflCoeff
    
    properties (Access=public)

    end
    
    methods
        function obj = PlaneWavePropagation(f)   % contructor
            obj@TransmAndReflCoeff(f); % invoking the superclass constructor
        end
        
        function p = planeWaveFluid(obj,zDistance)
            w = 2*pi*obj.f;
            u_z = obj.v0/i/w;
            p0 = -obj.rho_f*i^2*w^2*u_z/(-i*obj.h_f);    % Euler p. 120 KF
            p = p0.*exp(-i.*obj.h_f.*zDistance);
        end
    end
end



