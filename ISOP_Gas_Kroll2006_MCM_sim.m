function[S]=ISOP_Gas_Kroll2006_MCM_sim(test_dat,Do_GP)


% cd('/Users/cyshen/Documents/ucr/isoprene/F0AM-WAM-2020-distribute/Setups/Isoprene mechanism');
SZA=50;
expt=test_dat(1);
iC5H8 = test_dat(2);%ppb
    
iH2O2 = test_dat(3);
iH2O2 = iH2O2.*1000; %convert to ppb
%     
iNO = test_dat(4);%ppb
iNO2 = (test_dat(5)-iNO); %ppb

KT = test_dat(end)+273.13;
KP = 1013; %ASSUME 1ATM pressure

if Do_GP
    savename=['.\output\GP_Kroll2006Expt',...
        num2str(expt),'_MCM.mat'];
else
    savename=['.\output\Gas_Kroll2006Expt',...
        num2str(expt),'_MCM.mat'];
end
%% METEOROLOGY

Met = {...
%   names       values          
    'P'         KP; %Pressure, mbar
    'T'         KT; %Temperature, K
    'RH'        50; %Relative Humidity, percent
    'SZA'       SZA;
    'Nc'        1e3;
%     'Rp'        Rp*10^-4;  %cm
    'Cstar_threshold' 100;...   %Cstar threshold to filter compounds that can partition. ug/m3
    'kwall_loss'    0 ;... % default value of 2.5e-5;
    'kwall_lossV'   0 ;...
%     'LFlux'     'ExampleLightFlux.txt'; %Text file for radiation spectrum
    
    'EF_C5H8'   0  ;...          %emission factor to supply isoprene, continuous flow
    'EF_H2O2'   0  ;...          %emission factor to supply H2O2, continuous flow
    'EF_NO'     0  ;...  
    'J3'        0.000175./60;... % H2O2-OH+OH; Kroll state 0.00029/60,have used 0.00025/60 NOTE: 3.7e-6 s^-1 for Krechmer
    'J4'        0.29./60;... %was 0.29/60     %jvalue in Kroll seems too high (based on O3 production) Kroll 2006 %Override calcuated J-value for NO2 + hv -->
    'J5'        2E-4;...     %Override calculated J-value for NO3+hv--> NO   original=.02485
    'J6'        4.6E-4;...    %Override calculated J-value for NO3+hv--> NO2+O3   original=0.1747
    };
%% CHEMICAL CONCENTRATIONS
Nseed=1e3;
Rp=25e-7;% unit of cm, corresponds to 10nm
SAi=Nseed.*4*pi*Rp.^2;
Vseed=SAi.*Rp/3; % cm^2/cc
MWoa_init = 175;
Coa_ugm3 = .002;
iCoa = ugm3_to_mr(Coa_ugm3,MWoa_init).*1e9; %ppb. assume average MW of OA constituents 175g/mol

InitConc = {...
%   names       conc(ppb)           HoldMe
    'C5H8'      iC5H8           0;
    'H2O2'      iH2O2           0;
    'NO'        iNO             0;
    'NO2'       iNO2            0;
    'O3'        0               0;
    'CO'        0               1;   % ?
    'ttlOA'     iCoa            0;
    'OAinit'    iCoa            1;
    };

seedValues = {...
    'Rp'         Rp                  0;...    %initial seed radius (cm)
    'SA'         SAi                   0;...    %initial seed surface area (cm2/cm3)
    'Vseed'      Vseed                 0;...    %initial seed volume (cm3/cm3)
    'Nseed'      Nseed                 0;...    %initial number of seed particles
    'MWoa_init'   MWoa_init             1;...
    };
%% PARTICLE INITIALIZATION
%     Nseed = 1e3; %number seed per cm3
%     Rp = 25e-7; %seed particle radius (cm) (assuming 50nm diameter seed particles)
%     SAi = Nseed.*4*pi*(Rp)^2; %initial seed surface area cm2/cm3
%     Vseed = SAi.*Rp./3; %initial seed volume cm3/cm3
%     Dg = 0.05; % gas-phase diffusivity, cm2/s
%     MWoa_init = 175; %guess at average OA molecular mass (gets updated in model)
%     Coa_ugm3 = .002; %initial OA (ug/m3). going down to .001 throws errors, this is the lowest we can go    
%     gamma = 1;
%     %Note, if you change gamma (e.g. mass accommodation coefficient) you should also change ModelOptions.GPupdateTime to 10
%     %for gama = 1, recommed use ModelOptions.GPupdateTime = 10, otherwise computation errors result from stiff equations.
    
%% CHEMISTRY
ChemFiles = {...
    'MCMv331_K(Met)';
    'MCMv331_J(Met,0)'; %Jmethod flag of 0 specifies using "MCM" J-value method.
    'MCMv331_Inorg_Isoprene';
%     'Isoprene_DHDHP_scheme_V3_Thornton_ACSESC_2020';
    };
ParticleChemFiles= {...
%     'Isoprene_DHDHP_scheme_Particle_Chem'
    'MCMv331_particle_Chem';
    };%
    
%% Vapor Pressures // Molecular Mass // SMILES
% To run gas-particle partitioning module, you will need to have a
% saturation vapor concentration ("Cstar") for each specie you want to
% participate in the gas-particle partitioning. You can set a threshold
% Cstar above which the species is excluded from gas-particle partitioning.
% This helps the model run faster by limiting the number of species requiring the more intensive
% calculations

%If this is the first time running the model, you will need to find and
% possibly create two files (1) a .mat file containing info on each species
% "...SpeciesInfo.mat" and (2) an excel file containing the vapor pressures for
% each species. For the latter you may need to run
% GetSMILESforVaporPressure.m first which will create a third file (SMILES.txt). These files are usually in
% ...\Tools\SMILES\.

% matfile='C:\Users\tiann\desktop\scy\F0AM-WAM-2020-distribute\Tools\SMILES\MCMv331SpeciesInfo.mat';
% Vpfile = 'C:\Users\tiann\desktop\scy\F0AM-WAM-2020-distribute\Tools\SMILES\SMILES_vapPress_Compernolle.xlsx'; %full path to the vapor pressure file

matfile='..\F0AM-WAM-2020-distribute\Tools\SMILES\MCMv331SpeciesInfo.mat';
Vpfile = '..\F0AM-WAM-2020-distribute\Tools\SMILES\SMILES_vapPress_Compernolle.xlsx'; %full path to the vapor pressure file
%% DILUTION CONCENTRATIONS
% We are not diluting the chamber air, so this input is irrelevant (but still necessary).
BkgdConc = {...
%   names           values
    'DEFAULT'       0;   %0 for all zeros, 1 to use InitConc
    };

%%
wall_params.wcstar_thresh = 1E6; %cstar threshold to invoke wall partitioning
%equivalent absorbing organic concentration of the wall material C_wall
%in ug/m3. 
%Use C_wall = [] if you want to use the Krechmer et al parameterization based on c*'s. See F0AM_WallPartitioning_Generator.m
%Otherwise the C_wall value you enter here will be used for all
%species.  
wall_params.C_wall = [10000]; %30000; %equivalent wall absorbing mass concentration ug/m3  
wall_params.kwall_transport = 2.5e-4; %per second timescale to mix vapors to/from wall surface  default value 1e-5;
  

GPstruct.seedValues = seedValues;
GPstruct.gama = 1;
GPstruct.wall_params = wall_params;
GPstruct.Vpfile = Vpfile;
GPstruct.matfile = matfile;

%% OPTIONS
ModelOptions.Verbose       = 1;
ModelOptions.EndPointsOnly = 0;
ModelOptions.LinkSteps     = 0;
ModelOptions.IntTime       = 3600*10; % in seconds
ModelOptions.SavePath      = savename;
ModelOptions.GPupdateTime  = 60;
ModelOptions.SA            = 0;
ModelOptions.Do_Wall_Partitioning = 1;
% ModelOptions.GoParallel    = 0;
%% MODEL RUN
% SAPRC07tic_UCRisopc_output = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
% Do GP;
if Do_GP
    S=F0AM_ModelCore_D_GP(GPstruct,Met,InitConc,ChemFiles,ParticleChemFiles,BkgdConc,ModelOptions);
else
% Do gas phase;
    ModelOptions = rmfield(ModelOptions,{'GPupdateTime','Do_Wall_Partitioning','SA'});
%     Met=rmfield(Met,{'Nc'});
    rmi = strmatch('kwall_loss',Met(:,1)); 
    rmj = strmatch('Cstar_threshold',Met(:,1)); 
    rmk = strmatch('Nc',Met(:,1));
    rm = union(rmi,rmj);
    rm = union(rm,rmk);
    Met = Met(~ismember(1:length(Met),rm),:);
    S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);
end

%% IMPORTING DATA




