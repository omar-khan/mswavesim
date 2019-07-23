classdef mParams
    properties
        Ms;
        GAMMA;
        F;                       % Frequency 2 pi F/(GAMMA*Ms)
        ALPHA;
        sigma;
        mu0;
        alphae;
        Ha;
        L;
        Ex_Eff;                 % Exchange Coefficients (used in Tensor)
    end
    methods
        function obj = mParams(v1,v11,v2,v3,v4,v5,v6)
            obj.Ms = v1;
            obj.GAMMA = 2.21e5;
            obj.F = 2*pi*v2*1e9/(obj.GAMMA*obj.Ms);
            obj.ALPHA = v3;
            obj.sigma = v4;
            obj.mu0 = 4e-7*pi;
            obj.alphae = obj.sigma*v5.Vol_F(3)^2/12*obj.GAMMA*obj.Ms*obj.mu0;
           % obj.Ha = 2e5/(8*pi)/obj.Ms;
            obj.Ha = v11/obj.Ms;
            obj.L = v6*sqrt(2*1.3e-11/(obj.mu0*obj.Ms^2));
            obj.Ex_Eff = [obj.L^2/(2*v5.Cell_F(1)^2) obj.L^2/(2*v5.Cell_F(2)^2) obj.L^2/(2*v5.Cell_F(3)^2)];
        end
        function print(obj)
            displayBlock('Material Parameters');
            fprintf('Ha:                             %.2g A/m (%g T)\n', obj.Ha, obj.mu0*obj.Ha);
            fprintf('Ms:                             %.2g A/m (%g T)\n', obj.Ms, obj.mu0*obj.Ms);
            fprintf('Gamma:                          %.3g\n', obj.GAMMA);
            fprintf('Damping Parameter (alpha):      %.3f\n', obj.ALPHA);
            fprintf('Permeability Constant (mu0):    %g (=%g pi)\n', obj.mu0, obj.mu0/pi);
            fprintf('Exchange Length (m):            %g m\n', obj.L);
            fprintf('Eddy Current Damping (alpha_e): %g\n', obj.alphae);
            fprintf('Conductivity (sigma):           %g 1/(ohm m)\n', obj.sigma);
            fprintf('Frequency (W):                  %f (%.3g GHz)\n', obj.F, obj.F*obj.GAMMA*obj.Ms/2/pi/1e9);
        end
    end
end