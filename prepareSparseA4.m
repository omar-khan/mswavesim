function A = prepareSparseA4(geo,mat,BAND,N,m0,Ho,P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z)
    indice=0;
    melem=8*geo.B(1)*geo.B(2)*geo.B(3)*BAND(1)*BAND(2)*BAND(3);
   
   % mex mex_prepareSparseIndices.c
    %[ind1v,ind2v,indausv] = 
   % mex_prepareSparseIndices(geometry.B(1),geometry.B(2),geometry.B(3),BAND(1),BAND(2),BAND(3));
   
    fprintf('NP\n');
    indausv(1:melem)=0;
    tic
    for kk = 1:geo.B(3)
        for jj = 1:geo.B(2)
            for ii = 1:geo.B(1)
                for nn = max(1,kk-BAND(3)):min(geo.B(3),kk+BAND(3))
                    for mm = max(1,jj-BAND(2)):min(geo.B(2),jj+BAND(2))
                        for ll = max(1,ii-BAND(1)):min(geo.B(1),ii+BAND(1))
                            indice=indice+1;
                            indaus=(mod((ii-ll),2*geo.B(1))+1)+(mod((jj-mm),2*geo.B(2)))*geo.ZP(1)+(mod((kk-nn),2*geo.B(3)))*geo.ZP(1)*geo.ZP(2);
                            indausv(indice)=indaus;
                        end
                    end
                end
            end
        end
    end
    indausv=indausv(1:indice);
    NP=spblkdiag([reshape(N.xx(indausv),1,1,indice) reshape(N.xy(indausv),1,1,indice) reshape(N.xz(indausv),1,1,indice);reshape(N.xy(indausv),1,1,indice) reshape(N.yy(indausv),1,1,indice) reshape(N.yz(indausv),1,1,indice);reshape(N.xz(indausv),1,1,indice) reshape(N.yz(indausv),1,1,indice) reshape(N.zz(indausv),1,1,indice)]);
    clearvars indausv;
    toc
    
    fprintf('ind1,2\n');
    tic
    ind1v(1:indice)=0;
    ind2v(1:indice)=0;
    indice=0;
    for kk = 1:geo.B(3)
        for jj = 1:geo.B(2)
            for ii = 1:geo.B(1)
                ind1=(kk-1)*geo.B(2)*geo.B(1)+(jj-1)*geo.B(1)+ii; 
                for nn = max(1,kk-BAND(3)):min(geo.B(3),kk+BAND(3))
                    for mm = max(1,jj-BAND(2)):min(geo.B(2),jj+BAND(2))
                        for ll = max(1,ii-BAND(1)):min(geo.B(1),ii+BAND(1))
                            indice=indice+1;
                            ind2=(nn-1)*geo.B(2)*geo.B(1)+(mm-1)*geo.B(1)+ll; 
                            ind1v(indice)=ind1;
                            ind2v(indice)=ind2;
                        end
                    end
                end
            end
        end
    end
    ind1v=ind1v(1:indice);
    ind2v=ind2v(1:indice);
    toc
  
    tic
    P_l2a_z=spblkdiag(reshape(P_l2a_z(:,ind2v),3,1,indice));
    P_l2a_x=spblkdiag(reshape(P_l2a_x(:,ind2v),3,1,indice));
    P_a2l_z=spblkdiag(reshape(P_a2l_z(ind1v,:)',3,1,indice))';
    P_a2l_x=spblkdiag(reshape(P_a2l_x(ind1v,:)',3,1,indice))';
    toc
 
    fprintf('bc_index\n');
    tic
    for kk = 1:geo.B(3)
        for jj = 1:geo.B(2)
            for ii = 1:geo.B(1)
                if((geo.B(3) > 1) && (kk==1 || kk== geo.B(3)))
                   indbc=(kk-1)*geo.B(2)*geo.B(1)+(jj-1)*geo.B(1)+ii;
                   %corrections to NP
                   NP(indbc,(indbc-1)*9+1)=NP(indbc,(indbc-1)*9+1)+mat.Ex_Eff(3);
                   NP(indbc,(indbc-1)*9+5)=NP(indbc,(indbc-1)*9+5)+mat.Ex_Eff(3);
                   NP(indbc,(indbc-1)*9+9)=NP(indbc,(indbc-1)*9+9)+mat.Ex_Eff(3);
                end
                if((geo.B(2) > 1) && (jj==1 || jj== geo.B(2)))
                   indbc=(kk-1)*geo.B(2)*geo.B(1)+(jj-1)*geo.B(1)+ii;
                   %corrections to NP
                   NP(indbc,(indbc-1)*9+1)=NP(indbc,(indbc-1)*9+1)+mat.Ex_Eff(2);
                   NP(indbc,(indbc-1)*9+5)=NP(indbc,(indbc-1)*9+5)+mat.Ex_Eff(2);
                   NP(indbc,(indbc-1)*9+9)=NP(indbc,(indbc-1)*9+9)+mat.Ex_Eff(2);

                end
                if((geo.B(1) > 1) && (ii==1 || ii == geo.B(1)))
                   indbc=(kk-1)*geo.B(2)*geo.B(1)+(jj-1)*geo.B(1)+ii;
                   %corrections to NP
                   NP(indbc,(indbc-1)*9+1)=NP(indbc,(indbc-1)*9+1)+mat.Ex_Eff(1);
                   NP(indbc,(indbc-1)*9+5)=NP(indbc,(indbc-1)*9+5)+mat.Ex_Eff(1);
                   NP(indbc,(indbc-1)*9+9)=NP(indbc,(indbc-1)*9+9)+mat.Ex_Eff(1);
                end
            end
        end
    end
    toc
    
    fprintf('prepare A components\n');
    tic
    NPzz=P_a2l_z*NP*P_l2a_z;
    Azz = sparse(ind1v,ind2v,-diag(NPzz),geo.B(1)*geo.B(2)*geo.B(3),geo.B(1)*geo.B(2)*geo.B(3));
    clearvars NPzz;
    
    NPxz=P_a2l_x*NP*P_l2a_z;
    Axz = sparse(ind1v,ind2v,-diag(NPxz),geo.B(1)*geo.B(2)*geo.B(3),geo.B(1)*geo.B(2)*geo.B(3)); 
    clearvars NPxz;
    
    NPzx=P_a2l_z*NP*P_l2a_x;
    Azx = sparse(ind1v,ind2v,-diag(NPzx),geo.B(1)*geo.B(2)*geo.B(3),geo.B(1)*geo.B(2)*geo.B(3));
    clearvars NPzx;
    
    NPxx=P_a2l_x*NP*P_l2a_x;
    Axx = sparse(ind1v,ind2v,-diag(NPxx),geo.B(1)*geo.B(2)*geo.B(3),geo.B(1)*geo.B(2)*geo.B(3));
    clearvars NPxx;
    toc
    
    clearvars ind1v;
    clearvars ind2v;
    
    Azz=Azz+sparse(diag(Ho)) + i*mat.F*mat.ALPHA*speye(geo.B(1)*geo.B(2)*geo.B(3));
    Azx=Azx+(-i*mat.F)*speye(geo.B(1)*geo.B(2)*geo.B(3));
    Axz=Axz+i*mat.F*speye(geo.B(1)*geo.B(2)*geo.B(3));
    Axx=Axx+sparse(diag(Ho)) + i*mat.F*mat.ALPHA*speye(geo.B(1)*geo.B(2)*geo.B(3));
    
    A=[Azz Azx; Axz Axx];
end