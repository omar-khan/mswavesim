classdef cwgParams
properties
    width; 
    space;
    start;
    LW;     % Discretized Width
    LS;     % Discretized Spacing
    START;  % Discretized Starting Point
    CWG1;
    CWG2;
    CWG3;
    Jsimm;
    Jasimm;
    J;
end
methods
    function obj = cwgParams(v1,v2,v3,geo,v5,v6)
        % Wave Guide Geometry
        obj.width = v1;
        obj.space = v2;
        obj.start = v3;
        obj.Jsimm = v5;
        obj.Jasimm = v6;

        obj.LW    = floor(obj.width/geo.Cell_C(1));
        obj.LS    = floor(obj.space/geo.Cell_C(1));
        obj.START = floor(obj.start/geo.Cell_C(1));

        obj.CWG1 = [obj.START obj.START+obj.LW];
        obj.CWG2 = [obj.START+obj.LW+obj.LS obj.START+obj.LS+2*obj.LW];
        obj.CWG3 = [obj.START+2*obj.LS+2*obj.LW obj.START+2*obj.LS+3*obj.LW];

        % Prepare Waveguide ---------------------
        obj.J = field();
        obj.J.x = zeros(geo.BS(1),geo.BS(2),geo.BS(3));
        obj.J.y = zeros(geo.BS(1),geo.BS(2),geo.BS(3));
        obj.J.z = zeros(geo.BS(1),geo.BS(2),geo.BS(3));

        % Current Density Along Y
        obj.J.y(obj.CWG1(1):obj.CWG1(2),     1:geo.BS(2),    1:geo.BS(3)) = -obj.Jsimm/2-obj.Jasimm; % GND
        obj.J.y(obj.CWG2(1):obj.CWG2(2),     1:geo.BS(2),    1:geo.BS(3)) = obj.Jsimm;   % Source 
        obj.J.y(obj.CWG3(1):obj.CWG3(2),     1:geo.BS(2),    1:geo.BS(3)) = -obj.Jsimm/2+obj.Jasimm; % GND

    end
    function print(obj,geometry,material)
        displayBlock('CWG Parameters');
        fprintf('Current Density (J_Symm):       %g 1/m (%g kA/mm^2)\n', obj.Jsimm, obj.Jsimm*material.Ms*1e-9);
        fprintf('Current Density (J_Asymm):      %g 1/m (%g kA/mm^2)\n', obj.Jasimm, obj.Jasimm*material.Ms*1e-9);
        fprintf('Line Width:                     %.1f um\n', obj.LW*geometry.Cell_F(1)*1e6);
        fprintf('Line Spacing:                   %.1f um\n', obj.LS*geometry.Cell_F(1)*1e6);
        fprintf('Distance from L. Edge:          %.1f um\n', obj.START*geometry.Cell_F(1)*1e6);
        fprintf('Covered Area by CWG:            %.1f um\n', (2*obj.LS+3*obj.LW)*geometry.Cell_F(1)*1e6);
        fprintf('Distance from Film M. Point:    %.1f um (M.Point: %.0f um)\n', geometry.Vol_F(1)*1e6/2 - (obj.START+2*obj.LS+3*obj.LW)*geometry.Cell_F(1)*1e6, geometry.Vol_F(1)*1e6/2);
        fprintf('Current in Wave Guide:          %.1f\n',   geometry.Cell_F(1)*geometry.Cell_F(3)*sum(obj.J.y(:))/geometry.B(2));
    end
    function H = calculateCWGHa(obj, geo, C)
        sJ = field(); % Calculate Ha
        sJ.y = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
        sJ.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3)) = obj.J.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
        sJ.y = fftn(sJ.y);
        test = field();
        test.y = sJ.y;
        H = field(); 
        H.x = ifftn( test.y .* C.xy);
        H.z = ifftn(-test.y .* C.yz);
        H.x = H.x(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));        % Applied Field
        H.z = H.z(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));        % Applied Field
    end
    function Zeq = sourceAnalysis(obj,geo,mat,C,J,Hd,Ha)
        %% Move this part to new function
%         sJ = field(); % Calculate Ha
%         sJ.y = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
%         sJ.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3)) = obj.J.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
%         sJ.y = fftn(sJ.y);
%         test = field();
%         test.y = sJ.y;
%         H = field(); 
%         H.x = ifftn( test.y .* C.xy);
%         H.z = ifftn(-test.y .* C.yz);
%         H.x = H.x(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));        % Applied Field
%         H.z = H.z(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));        % Applied Field
        %       H.z = sum(sum(H.z+Hd.z,3),2);
%         H_temp = H;
        %% Move End
        
        x      = (1:geo.BS(1))*geo.Cell_F(1)-geo.Cell_F(1)/2;
        x      = x-(x(obj.CWG1(1))+x(obj.CWG3(2)))/2;
        cwg1_t = [obj.CWG1(1) obj.CWG1(2) obj.CWG2(1) obj.CWG2(2)];
        cwg2_t = [obj.CWG2(1) obj.CWG2(2) obj.CWG3(1) obj.CWG3(2)];

        % Total Flux
        delta_c  = geo.Cell_C(3) * geo.Cell_C(2) * geo.Cell_C(1);
        region_1 = sum(H.z(cwg1_t(1):cwg1_t(3)));
        region_2 = sum(H.z(cwg1_t(2):cwg1_t(4)));
        PHI1     = delta_c * region_1 + delta_c * region_2;

        Hzx      = H.z.*x(:);
        Hzxt1    = H.z.*(x(:)+geo.Cell_C(1)*(obj.LS+obj.LW/2));
        Hzxt2    = H.z.*(x(:)-geo.Cell_C(1)*(obj.LS+obj.LW/2));
        
        PHI1     = obj.LW*geo.Cell_F(1)/2*PHI1;
        region_3 = sum(Hzxt1(cwg1_t(1):cwg1_t(2)));
        region_4 = sum(Hzx(cwg1_t(3):cwg1_t(4)));
        PHI1     = PHI1 + delta_c*region_3 - delta_c * region_4;

        temp     = 0.5*i*PHI1*mat.mu0*mat.Ms*mat.Ms*mat.F*mat.GAMMA;
        POWER1   = temp*conj(obj.Jsimm/2+obj.Jasimm)*mat.Ms;
        U21      = temp*2/(obj.width*abs(geo.Vol_C(3)));
        
        region_5 = sum(H.z(cwg2_t(1):cwg2_t(3)));
        region_6 = sum(H.z(cwg2_t(2):cwg2_t(4)));
        PHI2     = -delta_c*region_5 - delta_c*region_6;
        PHI2     = obj.LW*geo.Cell_F(1)/2*PHI2;

        region_7 = sum(Hzx(cwg2_t(1):cwg2_t(2)));
        region_8 = sum(Hzxt2(cwg2_t(3):cwg2_t(4)));
        PHI2     = PHI2-delta_c*region_7 + delta_c*region_8;
        
        temp     = 0.5*i*PHI2*mat.mu0*mat.Ms*mat.Ms*mat.F*mat.GAMMA;
        POWER2   = temp*conj(obj.Jsimm/2-obj.Jasimm)*mat.Ms;
        U23      = temp*2/(obj.width*abs(geo.Vol_C(3)));

        POWER    = POWER1+POWER2;

        Zeq      = 2*POWER/(abs(obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3))))^2;
        
        %Check <average dissipation>:
        temp = 0.5 * sum( J.x.*conj(J.x) + J.z.*conj(J.z) ) * mat.mu0;
        temp*mat.ALPHA*delta_c*mat.F^2*mat.Ms^2*mat.GAMMA*mat.Ms;
    end
    function [Zmag Spower] = getImpedance(obj,geo,mat,J,Ha) % Works with 1D
        delta_f     = geo.Cell_F(3) * geo.Cell_F(2) * geo.Cell_F(1);
        in_current  = obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3));
        Spower      = (sum(J.x(:).*conj(Ha.x(:))+J.z(:).*conj(Ha.z(:)))*delta_f*0.5*mat.mu0*(i*mat.F*mat.GAMMA*mat.Ms)*mat.Ms^2);
        Zmag        = Spower*2/in_current^2; % Units Ohms
    end
    function Zmag = getImpedanceAnalytical(obj,geo,mat,Jc,Ha) % Works with 1D
        delta_f     = geo.Cell_F(1) * geo.Cell_F(2) * geo.Vol_F(3);
        in_current  = obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3));
        Spower      = (sum(Jc.x(:).*conj(Ha.x(:,:,1))+Jc.z(:).*conj(Ha.z(:,:,1)))*delta_f*0.5*mat.mu0*(i*mat.F*mat.GAMMA*mat.Ms)*mat.Ms^2);
        Zmag        = Spower*2/in_current^2; % Units Ohms
    end
    function Zmag = getImpedance2(obj,geo,mat,J,Ha) % Works with 1D
        delta_f     = geo.Cell_F(3) * geo.Cell_F(2) * geo.Cell_F(1);
        in_current  = obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3))*abs(geo.Vol_C(2));
        Spower      = (sum(J.x(:).*conj(Ha.x(:))+J.z(:).*conj(Ha.z(:)))*delta_f*0.5*mat.mu0*(i*mat.F*mat.GAMMA*mat.Ms)*mat.Ms^2);
        Zmag        = Spower*2/in_current^2; % Units Ohms
    end
    function Zmag = getImpedanceAnalytical2(obj,geo,mat,Jc,Ha) % Works with 1D
        delta_f     = geo.Cell_F(1) * geo.Cell_F(2) * geo.Vol_F(3);
        in_current  = obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3))*abs(geo.Vol_C(2));
        Spower      = (sum(Jc.x(:).*conj(Ha.x(:,:,1))+Jc.z(:).*conj(Ha.z(:,:,1)))*delta_f*0.5*mat.mu0*(i*mat.F*mat.GAMMA*mat.Ms)*mat.Ms^2);
        Zmag        = Spower*2/in_current^2; % Units Ohms
    end
    function [Z_21 Voltage] = pickupLineVoltage(obj,geo,mat,H,Hd)
        PHI1        = H.z(obj.CWG1(2):obj.CWG2(1),1:geo.BS(2),round(geo.BS(3)/2));
        PHI1        = PHI1 + Hd.z(obj.CWG1(2):obj.CWG2(1),1:geo.BS(2),round(geo.BS(3)/2));
        PHI1        = sum(PHI1(:))*mat.Ms*mat.mu0*geo.Cell_F(1)*geo.Cell_F(2);

        PHI2        = H.z(obj.CWG2(2):obj.CWG3(1),1:geo.BS(2),round(geo.BS(3)/2));
        PHI2        = PHI2 + Hd.z(obj.CWG2(2):obj.CWG3(1),1:geo.BS(2),round(geo.BS(3)/2));
        PHI2        = sum(PHI2(:))*mat.Ms*mat.mu0*geo.Cell_F(1)*geo.Cell_F(2);
        
    %    Voltage     = mat.F*mat.GAMMA*mat.Ms*abs(PHI1 - PHI2)/2;
        Voltage     = i*mat.F*mat.GAMMA*mat.Ms*(PHI1 - PHI2)/2;
        in_current  = obj.Jsimm*mat.Ms*obj.width*abs(geo.Vol_C(3));
        Z_21        = Voltage/in_current;
    end
end
end
