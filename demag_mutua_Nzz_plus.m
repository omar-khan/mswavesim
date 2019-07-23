% About: Inputs are:    ZP - Local Padding
%                       C  - Cell Sizes
%                       Xo,Yo,Zo - Offset Distances

function N=demag_mutua_Nzz_plus(ZP,C,Xo,Yo,Zo)

% Extra Padding
if ZP(1)>1; Nx=ZP(1)+2; else Nx=ZP(1)+3; end
if ZP(2)>1; Ny=ZP(2)+2; else Ny=ZP(2)+3; end
if ZP(3)>1; Nz=ZP(3)+2; else Nz=ZP(3)+3; end

% Debug FFT Inverse
iff = 1;
dx=C(1);
dy=C(2);
dz=C(3);

% Get F2's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F2=zeros(Nx,Ny,Nz);
F2=mex_demag_mutua_zz(F2,Nx,Ny,Nz,C(1),C(2),C(3),Xo,Yo,Zo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=0:(Nz/2)-1
%     k_pos=k+1;
%     k_inv=Nz-k;
%     k_mul=Zo+k*dz;
%     k_min=Zo-(k+1)*dz;
% 	for j=0:(Ny/2)-1
%         j_pos=j+1;
%         j_inv=Ny-j;
%         j_mul=Yo+j*dy;
%         j_min=Yo-(j+1)*dy;
%         for ii=0:(Nx/2)-1
%             i_pos=ii+1;
%             i_inv=Nx-ii;
%             i_mul=Xo+ii*dx;
%             i_min=Xo-(ii+1)*dx;
%             F2(i_pos,j_pos,k_pos) = NewellF2(k_mul,i_mul,j_mul);
% 			F2(i_inv,j_pos,k_pos) = NewellF2(k_mul,i_min,j_mul);
%             F2(i_pos,j_inv,k_pos) = NewellF2(k_mul,i_mul,j_min);
%             F2(i_inv,j_inv,k_pos) = NewellF2(k_mul,i_min,j_min);
%             F2(i_pos,j_pos,k_inv) = NewellF2(k_min,i_mul,j_mul);
%             F2(i_inv,j_pos,k_inv) = NewellF2(k_min,i_min,j_mul);
%             F2(i_pos,j_inv,k_inv) = NewellF2(k_min,i_mul,j_min);
% 			F2(i_inv,j_inv,k_inv) = NewellF2(k_min,i_min,j_min);
% 		end
% 	end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add 2 Extra Points
% Nx=Nx+2;
% Ny=Ny+2;
% Nz=Nz+2;

% Kronecker Product
Wx=sin(pi/Nx*(0:Nx-1)).^2;
Wy=sin(pi/Ny*(0:Ny-1)).^2;
% Wz=sin(pi/Nz*(0:Nz-1)).^2;
Wz=1-exp(1i*2*pi/Nz*(0:Nz-1));
W=reshape(kron(kron(Wz,Wy),Wx),Nx,Ny,Nz);

% Final Term
N=-16/(4*pi*dx*dy*dz)*(W.*fftn(F2));
if iff==1
	N=real(ifftn(N));
end

% Reductions
N(:,:,(Nz)/2+1) =[];
N(:,:,Nz/2)     =[];
N(:,(Ny)/2+1,:) =[];
N(:,Ny/2,:)     =[];
N(Nx/2+1,:,:)   =[];
N(Nx/2,:,:)     =[];
if (Nx-3==1) N=N(1,:,:); end
if (Ny-3==1) N=N(:,1,:); end
if (Nz-3==1) N=N(:,:,1); end