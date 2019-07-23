clc; 
clear variables;

source('config.m');

%% 2.   Prepare Classes
%% --------------------------------------------
geometry = mGeometry(Vol_C,Vol_F,Cell_C,Cell_F,FS_D);
material = mParams(Ms,Ha,F,ALPHA,sigma,geometry,L);
cwg_s    = cwgParams(cwg_width, cwg_space, cwg_starting_pos, geometry, Jsimm, Jasimm);
cwg_p 	 = cwgParams(cwg_width, cwg_space, cwg_starting_pos+3*cwg_width+2*cwg_space+cwg_pickup_distance, geometry, Jsimm, Jasimm);
sim      = mSim(rel,band,loop);
plotting = plotSettings('k','b',cwg_s,10);      % Plot Settings
FLIST = [F:F_Int:F_Last];

[m0,Ho] = prepareGroundState(geometry,material,init_m);
[P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z] = buildP_a2l(m0);

C    = tensor(geometry,0);  % 0 = Current Tensor (with 0 Distance)
C_f  = tensor(geometry,1);  % 1 = Current Tensor (Film Layer)
N    = tensor(geometry,2);  % 2 = Self Demag Tensor (Magnetic Film)
N_f  = tensor(geometry,3);  % 3 = Mutual Demag Tensor (Thin-Film2CWG)
N_b  = tensor(geometry,4);  % 4 = Demag Tensor (Single Cell in Z)

%% 5.  Print Things --------------------------
geometry.print();                               % Print Geometry
material.print();                               % Print Material Params
cwg_s.print(geometry,material);                 % print CWG Params
cwg_p.print(geometry,material); 				% Print Pickup Line params

H = getHaxHaz(geometry,cwg_s,C_f);
%[sph,abcs]=plotSpectralAppliedField(plotting,geometry,H,2);

%% 7.  Start Solver ---------------------------
count = 1;  
FofK = field();
FofK.x = zeros(geometry.B(1),length(FLIST));
FofK.z = zeros(geometry.B(1),length(FLIST));
for TT=FLIST
    material = mParams(Ms,Ha,TT,ALPHA,sigma,geometry,L);   % Prepare Material
    BAND = adjustBand(TT,band,1);                       % For Critical Frequencies (Inputs: Band + critical size)

    A = prepareSparseA4(geometry,material,BAND,N,m0,Ho,P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z);       % Prepare Sparse Tensor

    Precond = prepareSparseA4(geometry,material,[1 1 1],N,m0,Ho,P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z);       % Prepare Sparse Tensor

    [J,Zm(count),Voltage(count)] = config_perugiaCR19(Precond, geometry, material, cwg_s, cwg_p, plotting, H, C, N, N_f, A, loop, rel, Ho, P_a2l_x,P_a2l_z,P_l2a_x,P_l2a_z); %oringinal Zm=Z21
    [Z(count),Spower(count)] = cwg_s.getImpedance(geometry,material,J,H); % computes the power and Z11

    %sp = plotSpectralField(plotting,geometry,J,3,1);
    %FofK.x(:,count) = sp.x;
    %FofK.z(:,count) = sp.z;
   
    fprintf('Voltage: %g V\n', Voltage(count));
    count = count+1;
end