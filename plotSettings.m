classdef plotSettings
    properties
        col_1;
        col_2;
        PSTART;
        PEND;
    end
    methods
        function obj = plotSettings(v1,v2,cwg,v4)
            obj.col_1 = v1;
            obj.col_2 = v2;
            obj.PSTART = cwg.START-v4; 
            obj.PEND = cwg.START+2*cwg.LS+3*cwg.LW+v4;
%             obj.PSTART = cwg.START;
%             obj.PEND = cwg.START+2*cwg.LS+3*cwg.LW;
        end
    end
end

% why v4 was subtracted and added in lines 12 & 13 not clear