classdef mSim
    properties
        rel;
        band;
        loop;
        fig_no;
    end
    methods
        function obj = mSim(rel, band, loop)
            obj.rel  = rel;
            obj.band = band;
            obj.loop = loop;
            obj.fig_no = 1;
        end
        function obj = adjustBand(TT,crit)
%             if (TT >= 2.9) && (TT < 4) 
%                 obj.band = [obj.band(1)*crit, band(2), band(3)*crit]; % Adjust BAND
%             else                       BAND = band; end
        end
    end
end      