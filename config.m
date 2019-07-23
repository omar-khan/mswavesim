%% 1.  Parameters -----------------------------

%% 1a.  Geometry Parameters
%% --------------------------------------------
%% Geometry is configured as two layers. The top layer is contains 
%% the embedding of the source and pickup waveguides. The bottom layer
%% contains the embedding of the thin-film. The volume of the both
%% layers is given as Vol_C and Vol_F. Discrete cells in both layers
%% is given as Cell_C and Cell_F. FS_D is the distance between both
%% layers.
%% --------------------------------------------
Vol_C  = [50e-6,        50e-6,      120e-9];    % Volume (Waveguide)
Vol_F  = [Vol_C(1),     Vol_C(2),    30e-9];    % Volume (Thin Film)
Cell_C = [0.2e-6,       1e-6,        30e-9];    % Cell Size (Waveguide)
Cell_F = [0.2e-6,       1e-6,        30e-9];    % Cell Size (Thin Film)
FS_D   = 50e-9;                                 % Distance (Waveguide Thin-Film)

%% 1b.  Material Parameters 
%% --------------------------------------------
Ms      = 0.87/(4e-7*pi);   % Saturation Point
Ha      = 12000;            % Applied Field Strength
ALPHA   = 0.01;             % Alpha (dispersion)
sigma   = 1e6;              % Conductivity
L       = 0;                % Exchange Length (0|1 Disable|Enable)
F       = 4.3;              % Starting Frequency
F_Int   = 1;                % Interval Frequency
F_Last  = 4.3;              % Ending Frequency 
init_m  = [0 1 0];          % Initial Magnetization Distribution

%% 1c.  Waveguide (CWG) Parameters
%% --------------------------------------------
Jsimm               = 1.25; %  (0.2)/(130e-9*2e-6)/Ms %1.25;   m (A/m^2 normalized by A/m)
Jasimm              = 0;
cwg_width           = 2.31e-6;
cwg_space           = 1.26e-6;
cwg_starting_pos    = 15e-6;
cwg_pickup_distance = 2.8e-6;

%% 1d.  Simulation Parameters
%% --------------------------------------------
rel     = 1;                % Relaxation Parameter
band    = [12, 10, 2];      % For Sparse Matrix
loop    = 1000;             % Iteration loop steps

