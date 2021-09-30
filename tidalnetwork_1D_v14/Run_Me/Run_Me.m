%this is a script to run the model. ini_dia_*version* is the initialisation
%script to define the boundary conditions, grids and parameters. Boyo is
%the main script of this model.you can use this script if you want to do 
%some scenarios. All the required scripts for this model are:
%
% - ini_dia_*version*.m: define the parameters, channel configurations and boundary conditions here
% - ini_dia_matrix_*version*.m : script to make all matrices needed
% - ini_dia_sedtan_simpang_*version*.m : define the initial sediment transport if you
% want to simulate sediment transport and bedUpdate and do hot run.
% - Aquae_*version*.m : hidrodynamic calculation
% - sedtan_*version*.m : sediment transport calculation
% - simpang_*version*.m: nodal point calculation
% - bedte_*version*.m: bedupdate calculation
% - TiMor_*version*.m: bed update for tidally-averaged morfac (optional)
% - Boyo.m: the main script
% - harmfit.m: tidal component analysis in the post processing, if
% necessary

%% spatial grid size check
warning off
clear all
close all
clc
% BE AWARE THAT THE PATH IS NOT CORRECT FOR DIFFERENT VERSION!!
addpath(genpath('D:\tidalnetwork_1D_v14\source'));

ini_dia_3
Boyo_V8