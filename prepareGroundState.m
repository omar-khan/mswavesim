function [m0,Ho] = prepareGroundState(geometry, material, INIT)
    m0=zeros(geometry.B(1)*geometry.B(2)*geometry.B(3),3);
    m0(:,1)=INIT(1);
    m0(:,2)=INIT(2);
    m0(:,3)=INIT(3);
    
    Ho=zeros(geometry.B(1)*geometry.B(2)*geometry.B(3),3);
%     if INIT(1) > 0 Ho(:,1) = material.Ha; end;
%     if INIT(2) > 0 Ho(:,2) = material.Ha; end;
%     if INIT(3) > 0 Ho(:,3) = material.Ha; end;
    Ho=material.Ha*m0;
    
    Ho=dot(Ho',m0');
end