classdef tensor
    properties
        xx;
        xy;
        xz;
        yx;
        yy;
        yz;
        zx;
        zy;
        zz;
    end
    methods
        function obj = tensor(g,number)
            if (number == 0)
                mex mex_current_mutua_xy.cc;
                mex mex_current_mutua_yz.cc;
                mex mex_current_mutua_zx.cc;
                obj.xy = fftn(current_mutua_Cxy(g.ZPS,g.Cell_C,0,0,0)); % CS
                obj.yz = fftn(current_mutua_Cyz(g.ZPS,g.Cell_C,0,0,0)); % CS
            end
            if (number == 1)
                obj.xy = fftn(current_mutua_Cxy(g.ZPS,g.Cell_C,g.DS(1),g.DS(2),g.DS(3)));
                obj.yz = fftn(current_mutua_Cyz(g.ZPS,g.Cell_C,g.DS(1),g.DS(2),g.DS(3)));
                obj.zx = fftn(current_mutua_Czx(g.ZPS,g.Cell_C,g.DS(1),g.DS(2),g.DS(3)));
            end
            if (number == 2)
                mex mex_demag_mutua_xx.cc;
                mex mex_demag_mutua_xy.cc;
                mex mex_demag_mutua_xz.cc;
                mex mex_demag_mutua_yy.cc;
                mex mex_demag_mutua_yz.cc;
                mex mex_demag_mutua_zz.cc;
                obj.xx = (demag_mutua_Nxx(g.ZP,g.Cell_F,0,0,0));
                obj.xy = (demag_mutua_Nxy(g.ZP,g.Cell_F,0,0,0));
                obj.xz = (demag_mutua_Nxz(g.ZP,g.Cell_F,0,0,0));
                obj.yy = (demag_mutua_Nyy(g.ZP,g.Cell_F,0,0,0));
                obj.yz = (demag_mutua_Nyz(g.ZP,g.Cell_F,0,0,0));
                obj.zz = (demag_mutua_Nzz(g.ZP,g.Cell_F,0,0,0));
            end
            if (number == 3)
                obj.xz = fftn(demag_mutua_Nxz(g.ZPS,g.Cell_C,0,0,-g.DS(3))); 
                %sign of DS reversed because the CWG is now target
                obj.yz = fftn(demag_mutua_Nyz(g.ZPS,g.Cell_C,0,0,-g.DS(3)));
                obj.zz = fftn(demag_mutua_Nzz(g.ZPS,g.Cell_C,0,0,-g.DS(3)));
            end
            if (number == 4)
                ZP_a   = [g.ZP(1)     g.ZP(2)     1];
                Cell_a = [g.Cell_F(1) g.Cell_F(2) g.Vol_F(3)];
                obj.xx = fftn(demag_mutua_Nxx(ZP_a,Cell_a,0,0,0));
                obj.zz = fftn(demag_mutua_Nzz(ZP_a,Cell_a,0,0,0));
                obj.xz = fftn(demag_mutua_Nxz(ZP_a,Cell_a,0,0,0));
                obj.zx = obj.xz;
            end
            if (number == 5)
                mex mex_demag_mutua_xx.cc;
                mex mex_demag_mutua_xy.cc;
                mex mex_demag_mutua_xz.cc;
                mex mex_demag_mutua_yy.cc;
                mex mex_demag_mutua_yz.cc;
                mex mex_demag_mutua_zz.cc;
                obj.xx = (demag_mutua_Nxx_broken(g.ZP,g.Cell_F,0,0,0));
                obj.xy = (demag_mutua_Nxy_broken(g.ZP,g.Cell_F,0,0,0));
                obj.xz = (demag_mutua_Nxz_broken(g.ZP,g.Cell_F,0,0,0));
                obj.yy = (demag_mutua_Nyy_broken(g.ZP,g.Cell_F,0,0,0));
                obj.yz = (demag_mutua_Nyz_broken(g.ZP,g.Cell_F,0,0,0));
                obj.zz = (demag_mutua_Nzz_broken(g.ZP,g.Cell_F,0,0,0));            
            end
        end
        function obj = doFFT(v1)
            obj.xx = fftn(v1.xx);
            obj.xy = fftn(v1.xy);
            obj.xz = fftn(v1.xz);
            obj.yy = fftn(v1.yy);
            obj.yz = fftn(v1.yz);
            obj.zz = fftn(v1.zz);
        end
    end
end      