%% Parameterisation

%parameters and HD boundary conditions

name = 'tidal variation stability';
dt=300;                                       % time step in seconds.
store = 24*3600;                               %store the result
if store<dt
    ['warning you storing time step is smaller than the time step, so store = dt']
    ['lo nyimpen kebanyakan']
    store = dt;
elseif mod(store,dt)~=0
    ['warning the remainder after division of "store" by "dt" is not zero, change it so it is zero!']
    
end

dx= 250;                                % spatial step in meters; optimum dx = 2*B
Cd=2.725e-3;                            % Drag coefficient Chezy = 13sqrt(g) ( 60 =(2.725e-3))
g =9.81;                                % gravity
theta= 0.56;                            % time weight parameter for preissmann; >0.5 uncnditionally stable, the best accuracy between 0.55-0.6

%switch control
calc_Lb = 1;                % if efolding length scale needed = 1; no = 0
adv=1;                      %1 = include advection in HD calculation; 0 = no
calc_Q = 1;                 %1 = imposing equilibrium discharge; 0 self-defined discharge 
harm_anal = 1;              %1 = do harmonic analysis; 0 = no
harm_ini_final = 1;         %1 = do harmonic analysis for final bathymetry; 0 = no
calc_sedtan = 1;            %0 = no sediment calculation; 1 = yes
calc_bedte  = 1;            %0 = no bed update; 1 = yes
calc_TiMor  = 0;            %0 = normal morfac calculation; 1 = morfac for tide average sediment transport
switchtransp = 2;           % select sediment transport formulae, 1 = EH; 2 = vRijn 1984
switchndpt = 1;             % 0 = bolla, 1 = Arya, 2 = Wang
switchCp = 2;               % 0 = reference height depends on D50 (according to Einstein); 1 = reference height constant (1% of height) (minimum according to van Rijn 1984b); 2 = reference height depend on Cd, ks and depth (according to van Rijn 1984b); 
switchnpShi = 0;            % 0 = grain related Shields for nodal point relation, 1 = total Shields for nodal point, 2 = constant sqrt(Shields) (according to Baar et al 2018)
switchwetper = 1;           % 0 = changing wetted perimeter, 1 = constant wetted perimeter
switchstart = 0;            % 0 = cold start, 1 = hot start

if switchstart == 1
    coldmorf = load(''); % for bed level results from cold start simulation
    coldHD = load(''); % for HD and sediment transport results from cold start simulation
end
 
%sediment parameter
if switchtransp==1
    al_EH = 0.05;
    nt_EH = 2.5;
elseif switchtransp==2
    alb1 = 0.1;
    alb2 = 0.053;
    nt1 = 1.5;
    nt2 = 2.1;
    alI = 1;
end
rho_s = 2650;   % density of sediment
rho = 1000;     % density of water
Rr = (rho_s/rho)-1;     % spesific density of sediment
vk = 0.4;   % von karman
alc = (1/log10(exp(1))/vk);
ks = 0.035;           %roughness height (m)
tau_cr = 0.047;      %critical shear stress
temp = 15;          %°C
nu = 4e-5/(20+temp); %kinematic viscosity of water (1.2e-6 for cold water)
If = 1;				%Intermittency
lamp = 0.35;        %Bed Porosity
alw = 3;			%3, alpha parameter in Bolla NODAL POINT RELATION
rbolla = 0.5;        % calibrated parameter alphabn in ikeda 1982
morfac = 600;
morstt = 5*24*3600;   %bedupdate start

% parameters for nodal point relation methods
if switchndpt == 0 || switchndpt == 1
    alw = 3;			%3, alpha parameter in Bolla NODAL POINT RELATION
    rbolla = 0.5;        % calibrated parameter alphabn in ikeda 1982
elseif switchndpt == 2
    k_wang = 1;
end


%% INPUT sizes and topology - here you define the channel properties
% to add channel(s): just add line(s) of the matrices as the number of the channel that you want to add 

%Topo = (upstr1,upstr2,downstr1,downstr2) this is to define the connection of a certain channel with other channels
%upstream connections: (none if upstream boundary, one if bifur, two if confluence)
%downstream connections: (none if sea, one if confluence, two if bifurcation)
Topo = [...
    NaN NaN 2   3   ; ...1 -> indicating channel identity
    1   NaN NaN NaN ; ...2
    1   NaN NaN NaN ; ...3
    ];

%amount of channel(s) and boundaries - do not need to change 
Nc = length(Topo(:,1));
Nb = 4;  %number of boundary conditions prescribed for each channel (two conditions upstream (Q & Z) and two conditions downstream)    

%Sizes = [B0( width of the downstream end of each branch, if efolding length scale is considered otherwise the entire channel),
% L(channel length),Lb(efolding length scale),D50] for each branch

Sizes = [...
    1000 80000 1.5e4 0.00025; ...1
    500 20000 1.5e4 0.00025; ...2
    500 20000 1.5e4 0.00025; ...3    
    ];

Lc_up       = 40000         % the distance of channel widening upstream (from bifurcation) in meter


% In case the channel laength cannot be divided equally over dx
% this automatically modifies the channel length to discritize the length
% equally. do not need to change
for channel = 1:Nc
    if mod(Sizes(channel,2),dx)~=0
        test_length = round(Sizes(channel,2)/dx);
        Sizes(channel,2) = dx*test_length;
    end
end

%upstream end and downstream end of a channel (now all channel has the same
%slope (3e-5)
slope = [3e-5;...    1
         3e-5;...    2
         3e-5;...    3
         ];...    

% elevation of bed with respect to sea
eta_end = [Sizes(2,2)*slope(2);...    1
           0;...    2
           0;...    3
            ];...
            
 %slope correction for channel(s) that connect to other channels both downstream and usptream (if necessary)
%  slope = [2e-5;...      1
%          2e-5;...       2
%          2e-5;...       3
%          ];...   
 
% channel depth
H0=[10,...          1
    10,...         2
    10;...        3                
    ];
perturbdepth = 0;            %automatically change depth in both downstream channels...
                                %(one is deeper, the other shallower with this magnitude for hotstart run) 
%% Boundaries
% upstream
if calc_Q ==0
discharge = [0,...1
            ]; % Constant river discharge at landward boundary. self-defined one
end

%downstream
%tidal amplitude at system mouths
Tide_amp=[1.5 0;...    2
          1.5 0; ...   3
          ]; 
      
%  Tide_amp=[0 0;...    2
%           0 0; ...   3
%           ]; 
%tidal phase at system mouths in radian      
Tide_phase=[0 deg2rad(90);...    2
            0 deg2rad(90); ...   3
          ];

%tidal period in seconds
T_tide =[(12*3600);...
        (12*3600)/2;...
        ];      
jTide = T_tide(1)/dt; 
%% time discretation   
t_year = 6;                 %year
% t_day = 180;                %day
t_end = t_year*365*24*3600;  %duration, real one
% t_end = t_day*24*3600;  %duration, mainly for test simulation or short simulation
t=0:dt:t_end;           % time in seconds
Nt=length(t);

%% harmonic analysis in the middle of simulation
if harm_ini_final == 1
    harm_end  = [t_end];
    harm_step = 600;
    
    if harm_step<dt
        ['warning you harmonic analysis time step is smaller than the time step, so harm_step = dt']
        ['lo nyimpen kebanyakan']
    harm_step = dt;
    elseif mod(harm_step,dt)~=0
        ['warning the remainder after division of "harm_step" by "dt" is not zero, change it so it is zero!']
    
    end
end
