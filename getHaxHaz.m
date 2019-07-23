function [T] = getHaxHaz(geo,cwg_s,C)
    % Calculations:
    T = field();
    T.x = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    T.y = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    T.z = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    
    T.x(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3)) = cwg_s.J.x;
    T.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3)) = cwg_s.J.y;
    T.z(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3)) = cwg_s.J.z;

    T.x = fftn(T.x);
    T.y = fftn(T.y);
    T.z = fftn(T.z);
    
    test = field();
    test.x = T.x;
    test.y = T.y;
    test.z = T.z;
        
    T.x = test.y.*C.xy - test.z.*C.zx; 
    T.z = test.x.*C.zx - test.y.*C.yz;

    T.x=ifftn(T.x);
    T.z=ifftn(T.z);
   
    T.x=T.x(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
    T.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3))=0;
    T.z=T.z(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
    
    T.x=T.x(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3));
    T.y=T.y(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3));
    T.z=T.z(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3));
end