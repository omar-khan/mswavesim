function [sph,abcs] = plotSpectralAppliedField(plo,geo,H,number)
    sphx = field();
    sph.x=zeros(geo.ZP(1),geo.ZP(2),geo.ZP(3));
    sph.z=zeros(geo.ZP(1),geo.ZP(2),geo.ZP(3));
    sph.x(1:geo.B(1),1:geo.B(2),1:geo.B(3)) = H.x;
    sph.z(1:geo.B(1),1:geo.B(2),1:geo.B(3)) = H.z;
    sph.x=fftn(sph.x)/geo.B(1);
    sph.z=fftn(sph.z)/geo.B(1);
    h = figure(number)
    abcs = (0:geo.B(1)-1)*pi/geo.Vol_F(1)/1e6;
    plot(abcs, abs(sph.x(1:geo.B(1),1,1)),plo.col_1);
    hold on
    plot(abcs, abs(sph.z(1:geo.B(1),1,1)),plo.col_2);
    hold off
    legend('Hx','Hz');
    title('Spectral plot H_a','FontSize',16);
    xlabel('rad / cm 10^4')
    saveas(h,'plotSpectralAppliedField.eps','epsc');
end
