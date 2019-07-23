classdef mGeometry
    properties
        Vol_C;
        Vol_F;
        Cell_C;
        Cell_F;
        distanceFS;
        B;  % Number of Cells (Film)
        BS; % Number of Cells (CWG)
        DS;
        D;
        ZP;
        ZPS;
    end
    methods
        function obj = mGeometry(v1,v2,v3,v4,v5)
            obj.Vol_C = v1;
            obj.Vol_F = v2;
            obj.Cell_C = v3;
            obj.Cell_F = v4;
            obj.distanceFS = v5;

            obj.DS = [0, 0, -obj.distanceFS - obj.Vol_C(3)]; % Dot Displacement
            obj.D  = [obj.DS(1), obj.DS(2), -obj.distanceFS-obj.Vol_F(3)-obj.Vol_C(3)/2+obj.Cell_C(3)];    % Dot Displacement

            % Number of cells
            obj.B  = [round(obj.Vol_F(1)/obj.Cell_F(1)), round(obj.Vol_F(2)/obj.Cell_F(2)), round(obj.Vol_F(3)/obj.Cell_F(3))];
            obj.BS = [round(obj.Vol_C(1)/obj.Cell_C(1)), round(obj.Vol_C(2)/obj.Cell_C(2)), round(obj.Vol_C(3)/obj.Cell_C(3))];

            % How much local padding? 
            if obj.B(1)==1 obj.ZP(1)=obj.B(1); else obj.ZP(1)=round(obj.B(1)*2); end
            if obj.B(2)==1 obj.ZP(2)=obj.B(2); else obj.ZP(2)=round(obj.B(2)*2); end
            if obj.B(3)==1 obj.ZP(3)=obj.B(3); else obj.ZP(3)=round(obj.B(3)*2); end

            % How much local padding? (for Current)
            if obj.BS(1)==1 obj.ZPS(1)=obj.BS(1); else obj.ZPS(1)=round(obj.BS(1)*2); end
            if obj.BS(2)==1 obj.ZPS(2)=obj.BS(2); else obj.ZPS(2)=round(obj.BS(2)*2); end
            if obj.BS(3)==1 obj.ZPS(3)=obj.BS(3); else obj.ZPS(3)=round(obj.BS(3)*2); end
        end
        
        function print(obj)
            displayBlock('Geometry Parameters');
            fprintf('Total Volume (Source):          (%d x %d) um x %g nm\n', obj.Vol_C(1)*1e6,obj.Vol_C(2)*1e6,obj.Vol_C(3)*1e9);
            fprintf('Total Volume:                   (%d x %d) um x %g nm\n', obj.Vol_F(1)*1e6,obj.Vol_F(2)*1e6,obj.Vol_F(3)*1e9);
            fprintf('Discretized Cell (Source):      (%d x %d x %d)\n', obj.BS(1),obj.BS(2),obj.BS(3));
            fprintf('Discretized Cells:              (%d x %d x %d)\n', obj.B(1),obj.B(2),obj.B(3));
            fprintf('Flux Plane-Thin Film distance:  %g nm\n', abs(obj.D(3)*1e9));
            fprintf('Displacement:                   %g nm\n', abs(obj.DS(3)*1e9));
        end
    end
end      
