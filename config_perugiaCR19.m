function [Ja Z_21 Voltage] = config_perugiaCR19(Precond,geo, mat, cwg_s, cwg_p, plotting, H, C, N, N_f, A, loop, rel, Heff, P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z)
    %% 2. ...
    Ha = H;

    H.y((1:geo.B(1)),(1:geo.B(2)),(1:geo.B(3)))=0;  %  to be changed in future!!
    Hxla=H.x.*reshape(P_a2l_x(:,1),geo.B(1),geo.B(2),geo.B(3))+H.y.*reshape(P_a2l_x(:,2),geo.B(1),geo.B(2),geo.B(3))+H.z.*reshape(P_a2l_x(:,3),geo.B(1),geo.B(2),geo.B(3));
    Hzla=H.x.*reshape(P_a2l_z(:,1),geo.B(1),geo.B(2),geo.B(3))+H.y.*reshape(P_a2l_z(:,2),geo.B(1),geo.B(2),geo.B(3))+H.z.*reshape(P_a2l_z(:,3),geo.B(1),geo.B(2),geo.B(3));

   rhs=[Hzla(:); Hxla(:)];
     rhs=A\rhs;
 %   [rhs,x]=gmres(A,rhs,10,1e-6,10,Precond);
  %  if (x ~= 0) fprintf('gmres error'); end
    J.z=rhs(1:geo.B(1)*geo.B(2)*geo.B(3));
    J.x=rhs(geo.B(1)*geo.B(2)*geo.B(3)+1:2*geo.B(1)*geo.B(2)*geo.B(3));
    
    %% 3. Append Dispersion Matrix to Tensor ----
    % N = appendDispersion(N,mat);
        N.xx = fftn(N.xx);
        N.xy = fftn(N.xy);
        N.xz = fftn(N.xz);
        N.yy = fftn(N.yy);
        N.yz = fftn(N.yz);
        N.zz = fftn(N.zz);

    %% 4. Run Iterative Solver ------------------
    displayBlock('Error (E)');
    Ja = field();
    for k = 1:loop
        J.x=reshape(J.x,geo.B(1),geo.B(2),geo.B(3));
        J.z=reshape(J.z,geo.B(1),geo.B(2),geo.B(3));
        Jold=J;
        
        % L2A
        Ja.x  = P_a2l_x(:,1) .* J.x(:) + P_a2l_z(:,1) .* J.z(:);
        Ja.y  = P_a2l_x(:,2) .* J.x(:) + P_a2l_z(:,2) .* J.z(:);
        Ja.z  = P_a2l_x(:,3) .* J.x(:) + P_a2l_z(:,3) .* J.z(:);
        
        Ja.x=reshape(Ja.x,geo.B(1),geo.B(2),geo.B(3));
        Ja.y=reshape(Ja.y,geo.B(1),geo.B(2),geo.B(3));
        Ja.z=reshape(Ja.z,geo.B(1),geo.B(2),geo.B(3));
    
        Ha1 = getHd(geo,N,Ja);
        Ha1 = exchg_bc(Ha1,Ja,geo,mat);

        % A2L
        H1 = field();
        H1.x=Ha1.x.*reshape(P_a2l_x(:,1),geo.B(1),geo.B(2),geo.B(3))+Ha1.y.*reshape(P_a2l_x(:,2),geo.B(1),geo.B(2),geo.B(3))+Ha1.z.*reshape(P_a2l_x(:,3),geo.B(1),geo.B(2),geo.B(3));
        H1.z=Ha1.x.*reshape(P_a2l_z(:,1),geo.B(1),geo.B(2),geo.B(3))+Ha1.y.*reshape(P_a2l_z(:,2),geo.B(1),geo.B(2),geo.B(3))+Ha1.z.*reshape(P_a2l_z(:,3),geo.B(1),geo.B(2),geo.B(3));

        % Get Hxla
        E.x=Hxla(:)+H1.x(:)-(Heff(:)+i*mat.F*mat.ALPHA).*J.x(:)-i*mat.F*J.z(:);
        E.z=Hzla(:)+H1.z(:)-(Heff(:)+i*mat.F*mat.ALPHA).*J.z(:)+i*mat.F*J.x(:);
%         E.x=Hxla(:)+H1.x(:);
%         E.z=Hzla(:)+H1.z(:);
        Err(k)=max(abs(E.x(:)).^2+abs(E.z(:)).^2);
        
        rhs=[E.z(:); E.x(:)];
        clearvars E;
         rhs=A\rhs;
%        [rhs,x]=gmres(A,rhs,10,1e-6,10,Precond);
%        if (x ~= 0) fprintf('gmres error'); end

        
        J.z=Jold.z+rel*reshape(rhs(1:geo.B(1)*geo.B(2)*geo.B(3)),geo.B(1),geo.B(2),geo.B(3));
        J.x=Jold.x+rel*reshape(rhs(geo.B(1)*geo.B(2)*geo.B(3)+1:2*geo.B(1)*geo.B(2)*geo.B(3)),geo.B(1),geo.B(2),geo.B(3));
        if (k > 5)
           if ((Err(k) > Err(k-4)) || (Err(k)<1e-40))
               break;
           end
        end
        fprintf('%d\t', k);
        fprintf('%g\n', Err(k));
    end
   
    Ja.x  = P_a2l_x(:,1) .* J.x(:) + P_a2l_z(:,1) .* J.z(:);
    Ja.y  = P_a2l_x(:,2) .* J.x(:) + P_a2l_z(:,2) .* J.z(:);
    Ja.z  = P_a2l_x(:,3) .* J.x(:) + P_a2l_z(:,3) .* J.z(:);
    
    Ja.x=reshape(Ja.x,geo.B(1),geo.B(2),geo.B(3));
    Ja.y=reshape(Ja.y,geo.B(1),geo.B(2),geo.B(3));
    Ja.z=reshape(Ja.z,geo.B(1),geo.B(2),geo.B(3));
    
    %% 5. Calculating Mutual Demag Field (Hd) ----------
    H.x = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    H.y = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    H.z = zeros(geo.ZPS(1),geo.ZPS(2),geo.ZPS(3));
    H.x(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3)) = Ja.x;
    H.z(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3)) = Ja.z;
    H.y(1:geo.BS(1),1:geo.BS(2),geo.BS(3)-geo.B(3)+1:geo.BS(3)) = Ja.y;

    
    H.x=fftn(H.x);
    H.z=fftn(H.z);
    H.y=fftn(H.y);

    test = field();
    test.x=H.x;
    test.z=H.z;
    test.y=H.y;

    H.x =  test.x.*N_f.xz;
    H.z =  test.z.*N_f.zz;
    H.y =  test.z.*N_f.yz;

    H.x = ifftn(H.x);
    H.z = ifftn(H.z);
    H.y = ifftn(H.y);

    Hd = field();
    Hd.x = H.x(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));        % Magnetostatic Field
    Hd.z = H.z(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
    Hd.y = H.y(1:geo.BS(1),1:geo.BS(2),1:geo.BS(3));
    
    Hd.z = Hd.z+Hd.x+Hd.y;

    Hc = cwg_s.calculateCWGHa(geo, C);
    
    [Z_21 Voltage] = cwg_p.pickupLineVoltage(geo, mat, Hc, Hd);
end
