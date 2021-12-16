
%%                                                  MANTLE DIAGNOSTICS 2.33
%
% calculates connected areas by accounting for MD.addFactor times the area of
% the original plot: one part added at x=0 and one other part added at x=end
%
% evaluates the residual temperature
% field indicating +1 for hot plumes and -1 for cold plumes and 0 for the residual field
% ( hot plumes are defined by e.g. Tmean+0.5*(Tmax-Tmean) )
%
%                                                Fabio Crameri, 23.03.2020
% 						 Anna Guelcher, 10.10.2021
% 						 --> diagnostics as explained in 2022 G^3 paper on
% 						     strain-weakening rehology in Earth's lower mantle

function [MANTLE] = f_DiagnosticsMantle(FILE,GRID,SETUP,SWITCH,PLOT,STYLE)
%% DEFAULTS
MD.dummy = [];
if ~isfield(MD,'plumeDefinition');            MD.plumeDefinition            = 1;            end     %1: Solely by temperature, 2: Temperature & radial velocity
if ~isfield(MD,'activeUpwellingDefinition')   MD.activeUpwellingDefinition  = 2;            end     %1: Labrosse (2002,EPSL),  2: vz_res*T_res & T_res threshold !AG! 
if ~isfield(MD,'activeDownwellingDefinition') MD.activeDownwellingDefinition= 2;            end     %1: Labrosse (2002,EPSL),  2: vz_res*T_res & T_res threshold !AG!
if ~isfield(MD,'plotPLUMES');                 MD.plotPLUMES                 = logical(1); 	end     %plot comparison of all anomalies and selected plumes
if ~isfield(MD,'addFactor');                  MD.addFactor                  = 0.25;         end     %percent of total x-extend that is added at x=1 and x=end
if ~isfield(MD,'Tresidualmode');              MD.TresidualMode              = 1;            end     %1: horizontal, 2: horizontal band, 3: global, 4: regional
if ~isfield(MD,'DepthThresholdMode');         MD.DepthThresholdMode         = 2;            end     %1: Upper & lower z-threshold a plume/slab must both cross, 2: dz threshold !AG!
if ~isfield(MD,'upperDepthThresholdHot'); 	  MD.upperDepthThresholdHot 	= 0.7;          end     %upper threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'lowerDepthThresholdHot'); 	  MD.lowerDepthThresholdHot 	= 0.8;          end     %lower threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'upperDepthThresholdCold');	  MD.upperDepthThresholdCold	= 0.08;         end     %upper threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'lowerDepthThresholdCold');	  MD.lowerDepthThresholdCold	= 0.10;     	end     %lower threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'dzThresholdHot');             MD.dzThresholdHot             = 0.25;         end     %dz threshold PLUMES !AG!
if ~isfield(MD,'dzThresholdCold');            MD.dzThresholdCold            = 0.2;          end     %dz threshold SLABS  !AG!
if ~isfield(MD,'depthLevelDiagnostics');      MD.depthLevelsMode            = 1;            end     % same (1) or different (2)  UM, MM, LM depth definition for plumes/slabs !AG!

%% DEFAULT OUTPUT FIELDS
MANTLE.DiagnosticsFailed        = false;
MANTLE.upWelling                = NaN;      %1 for upwelling, 0 else (defined by threshold)
MANTLE.downWelling              = NaN;      %1 for downwelling, 0 else (defined by threshold)
MANTLE.upWellingAbsolute      	= NaN;      %1 for upwelling, 0 else
MANTLE.downWellingAbsolute   	= NaN;      %1 for downwelling, 0 else
MANTLE.plumesHot                = NaN;      %1 for hot plumes, 0 else
MANTLE.plumesCold               = NaN;  	%1 for cold plumes, 0 else
MANTLE.numHotPlumes             = NaN;    	%num
MANTLE.numColdPlumes            = NaN;  	%num
MANTLE.UpwellingVolume          = NaN;   	%[nd] or [m^2]or[m^3]
MANTLE.DownwellingVolume        = NaN;  	%[nd] or [m^2]or[m^3]
MANTLE.UpwellingVolPerc         = NaN;      %in [% of model domain]
MANTLE.DownwellingVolPerc       = NaN;  	%in [% of model domain]

MANTLE.VHmaxUM                 	= NaN;    	%[plotting dim]
MANTLE.VHmaxMM                 	= NaN;    	%[plotting dim]
MANTLE.VHmaxLM                 	= NaN;    	%[plotting dim]
MANTLE.VHmeanUM               	= NaN;    	%[plotting dim]
MANTLE.VHmeanMM                	= NaN;    	%[plotting dim]
MANTLE.VHmeanLM               	= NaN;    	%[plotting dim]

MANTLE.activeUpwelling          = NaN;  	%1 for active upwelling
MANTLE.activeDownwelling        = NaN;   	%1 for active downwelling
MANTLE.numActUpwelling          = NaN;    	%num
MANTLE.numActDownwelling        = NaN;   	%num
MANTLE.ActUpwellingVolume       = NaN;     	%[nd] or [m^2]or[m^3]
MANTLE.ActDownwellingVolume    	= NaN;     	%[nd] or [m^2]or[m^3]
MANTLE.ActUpwellingVolPerc      = NaN;   	%in [% of total upwelling]
MANTLE.ActDownwellingVolPerc  	= NaN;    	%in [% of total downwelling]

MANTLE.passiveUpwelling         = NaN;    	%1 for passive upwelling
MANTLE.passiveDownwelling       = NaN;   	%1 for passive downwelling
MANTLE.numPassUpwelling         = NaN;    	%num
MANTLE.numPassDownwelling       = NaN;    	%num
MANTLE.PassUpwellingVolume   	= NaN;    	%[nd] or [m^2]or[m^3]
MANTLE.PassDownwellingVolume 	= NaN;      %[nd] or [m^2]or[m^3]
MANTLE.PassUpwellingVolPerc  	= NaN;  	%in [% of total upwelling]
MANTLE.PassDownwellingVolPerc 	= NaN;      %in [% of total downwelling]

%!AG! output if diagnostics within volume domains
% active Upwellings               [Min, Max, Mean]
% if seperate for each plumes: 2nd dimension: MANTLE.numActUpwelling  
MANTLE.ActUpwellingT         = [NaN, NaN, NaN];
MANTLE.ActUpwellingTres      = [NaN, NaN, NaN];
MANTLE.ActUpwellingEta       = [NaN, NaN, NaN];
MANTLE.ActUpwellingVstr      = [NaN, NaN, NaN];
MANTLE.ActUpwellingVz        = [NaN, NaN, NaN];
% active downwellings             [Min, Max, Mean]
% if seperate for each plumes: 2nd dimension: MANTLE.numActDownwelling  
MANTLE.ActDownwellingT       = [NaN, NaN, NaN];
MANTLE.ActDownwellingTres    = [NaN, NaN, NaN];
MANTLE.ActDownwellingEta     = [NaN, NaN, NaN];
MANTLE.ActDownwellingVstr    = [NaN, NaN, NaN];
MANTLE.ActDownwellingVz      = [NaN, NaN, NaN];

MANTLE.plumeHotNumberUM         = NaN;    	%num
MANTLE.plumeHotVHmaxUM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminUM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanUM         = NaN;    	%[plotting dim]
MANTLE.plumeHotNumberMM         = NaN;    	%num
MANTLE.plumeHotVHmaxMM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminMM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanMM         = NaN;    	%[plotting dim]
MANTLE.plumeHotNumberLM         = NaN;    	%num
MANTLE.plumeHotVHmaxLM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminLM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanLM         = NaN;    	%[plotting dim]

MANTLE.plumeColdNumberUM      	= NaN;    	%num
MANTLE.plumeColdVHmaxUM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminUM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanUM    	= NaN;    	%[plotting dim]
MANTLE.plumeColdNumberMM    	= NaN;    	%num
MANTLE.plumeColdVHmaxMM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminMM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanMM        = NaN;    	%[plotting dim]
MANTLE.plumeColdNumberLM     	= NaN;    	%num
MANTLE.plumeColdVHmaxLM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminLM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanLM        = NaN;    	%[plotting dim]

MANTLE.contLocationX        	= NaN;    	%[plotting dim]
MANTLE.contLocationZ            = NaN;      %[plotting dim]

MANTLE.contNumberUM             = NaN;    	%num
MANTLE.contVHmaxUM              = NaN;    	%[plotting dim]
MANTLE.contVHminUM              = NaN;    	%[plotting dim]
MANTLE.contVHmeanUM             = NaN;    	%[plotting dim]


MANTLE.llsvp                    = NaN;          % 1 for llsvp !AG!	
MANTLE.llsvpVolume              = 0;            %[nd] or [m^2]or[m^3] !AG!
MANTLE.llsvpVolPerc		= 0;            %in [% of total mantle] !AG!
MANTLE.llsvpCMBPerc             = 0;            %in [% of CMB coverage] !AG!
MANTLE.llsvpLocationX        	= NaN;    	%[plotting dim]
MANTLE.llsvpLocationZ           = NaN;          %[plotting dim]

MANTLE.llsvpNumberLM            = NaN;    	%num
MANTLE.llsvpVHmaxLM             = NaN;    	%[plotting dim]
MANTLE.llsvpVHminLM             = NaN;    	%[plotting dim]
MANTLE.llsvpVHmeanLM            = NaN;    	%[plotting dim]

if strcmp(GRID.Type,'yinyang')
    %add here......
end

%% INPUT ADJUSTMENTS
if strcmp(GRID.Type,'yinyang') %grid might NOT be flipped in yinyang
    warning('check if these values need to be flipped upside down!');
%     MD.upperDepthThresholdHot  = 1-MD.upperDepthThresholdHot;
%     MD.lowerDepthThresholdHot  = 1-MD.lowerDepthThresholdHot;
%     MD.upperDepthThresholdCold = 1-MD.upperDepthThresholdCold;
%     MD.lowerDepthThresholdCold = 1-MD.lowerDepthThresholdCold;
end

%% SETUP SWITCHES
if strcmp(GRID.Type,'yinyang')
    nb = 2;
else
    nb = 1;
end

%% READ INPUT DATA
%temperature
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Temperature';
DATA.FieldAbbreviation      = 'T';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
interiorIndexLow            = floor( size(GRID.Z_3D,3)*1/10 );  %bottom of the mantle
interiorIndexHigh           = ceil( size(GRID.Z_3D,3)*8/10 );   %top of the mantle
if strcmp(GRID.Type,'yinyang')
    T_3D        = PLOT.T_3Dyin; T_3Dyang   = PLOT.T_3Dyang;
    Tmin        = min([PLOT.T_3Dyin(:);PLOT.T_3Dyang(:)]);
    Tmax        = max([PLOT.T_3Dyin(:);PLOT.T_3Dyang(:)]);
    dummy       = [PLOT.T_3Dyin(:,:,interiorIndexLow:interiorIndexHigh);PLOT.T_3Dyang(:,:,interiorIndexLow:interiorIndexHigh)];
    Tmedian    	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
else %all other grid types
    T_3D        = PLOT.T_3D;
    Tmin        = min(PLOT.T_3D(:));
    Tmax        = max(PLOT.T_3D(:));
    dummy       = PLOT.T_3D(:,:,interiorIndexLow:interiorIndexHigh);
    Tmedian  	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
end
clearvars dummy
if DATA.NotFound
    MANTLE.DiagnosticsFailed	= true;
    warning off backtrace
    disp(' '); warning('Mantle diagnostics failed due to missing temperature data.'); disp(' ');
    warning on backtrace
    return
    
end

%continental composition
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Cont. crust';
DATA.FieldAbbreviation      = 'CC';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    CC_3D    	= PLOT.CC_3Dyin; CC_3Dyang   = PLOT.CC_3Dyang;
else %all other grid types
    CC_3D  	= PLOT.CC_3D;
end
clearvars dummy

%primordial composition
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Primordial';
DATA.FieldAbbreviation      = 'PRM';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    PRM_3D    	= PLOT.PRM_3Dyin; PRM_3Dyang   = PLOT.PRM_3Dyang;
else %all other grid types
    PRM_3D  	= PLOT.PRM_3D;
end
clearvars dummy

%velocity
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Velocity';
DATA.FieldAbbreviation      = 'V';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
%     VX_3D = PLOT.VX_3Dyin; VX_3Dyang = PLOT.VX_3Dyang;
%     VY_3D = PLOT.VY_3Dyin; VY_3Dyang = PLOT.VY_3Dyang;
    VH_3D = PLOT.VH_3Dyin; VH_3Dyang = PLOT.VH_3Dyang;
    VZ_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
    VZmax       = max(abs([VZ_3D(:);VZ_3Dyang(:)]));   %doesn't take variable gridspacing into account....
    VZupmax     = abs(max([VZ_3D(:);VZ_3Dyang(:)]));
    VZdownmax 	= abs(min([VZ_3D(:);VZ_3Dyang(:)]));    
else %all other grid types
%     VX_3D       = PLOT.VX_3D;
%     VY_3D       = PLOT.VY_3D;
    VH_3D       = PLOT.VH_3D;
    VZ_3D       = PLOT.VZ_3D;
    VZmax       = max(abs(VZ_3D(:)));   %doesn't take variable gridspacing into account....
    VZupmax   	= abs(max(VZ_3D(:)));
    VZdownmax  	= abs(min(VZ_3D(:)));
end
if DATA.NotFound
    MANTLE.DiagnosticsFailed	= true;
    warning off backtrace
    disp(' '); warning('Mantle diagnostics failed due to missing velocity data.'); disp(' ');
    warning on backtrace
    return
end

%!AG! add basalt, viscosity and viscous strain field, for mantle domain diagnostics
% basalt
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Basalt';
DATA.FieldAbbreviation      = 'BS';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    BS_3D      = PLOT.BS_3Dyin; ETA_BSyang   = PLOT.BS_3Dyang;
else %all other grid types
    BS_3D      = PLOT.BS_3D;
end
clearvars dummy

% viscosity
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Viscosity';
DATA.FieldAbbreviation      = 'ETA';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    ETA_3D    	= PLOT.ETA_3Dyin; ETA_3Dyang   = PLOT.ETA_3Dyang;
else %all other grid types
    ETA_3D  	= PLOT.ETA_3D;
end
clearvars dummy

% viscous strain
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Viscous strain';
DATA.FieldAbbreviation      = 'VSTR';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    VSTR_3D    	= PLOT.VSTR_3Dyin; ETA_3Dyang   = PLOT.VSTR_3Dyang;
else %all other grid types
    VSTR_3D  	= PLOT.VSTR_3D;
end
if DATA.NotFound
    disp(' '); warning('No viscous strain field.'); disp(' ');
end
clearvars dummy


%% DIMENSIONALISATION
VH_3D           = VH_3D*SETUP.Vscale;
if strcmp(GRID.Type,'yinyang')
    VH_3Dyang 	= VH_3Dyang*SETUP.Vscale;
end

%% EXTRACT FLOW-VELOCITY CHARACTERISTICS AT CERTAIN DEPTH LEVELS
DomainDepthMax          = max(GRID.Z_3Dp(:));
DomainDepthMin          = min(GRID.Z_3Dp(GRID.Z_3Dp>0));
DomainDepth             = abs( DomainDepthMax - DomainDepthMin );
%define depth levels
if (MD.depthLevelsMode == 1)
    UMdepth       	= 0.85;               %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    MMdepth             = 0.48;                %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    LMdepth             = 0.22;               %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    %!AG! why is this flipped with respect to (1: bottom, 0: top)?
    DepthIdxHot  = 1;
    DepthIdxCold = 1;
elseif (MD.depthLevelsMode == 2) %!AG!
    UMdepth             = [0.85; 0.9];         %PLUME; SLAB, 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    MMdepth             = [0.5;  0.7];         %PLUME; SLAB, 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    LMdepth             = [0.2;  0.5];         %PLUME; SLAB, 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
    
    DepthIdxHot  = 1;
    DepthIdxCold = 2;
end

depthLevelDiagnostics       = [UMdepth, MMdepth, LMdepth];  %values between 1 (surface) and 0 (bottom)
%derive actual depth values and indices
depthLevelDiagnosticsP      = DomainDepthMax - depthLevelDiagnostics.*DomainDepth;  %plotting values
depthLevelDiagnosticsIdx    = zeros(size(depthLevelDiagnosticsP)); %array indices

for idD=1:size(depthLevelDiagnosticsP,2) % loop over UM, MM, LM 
    for idP=1:size(depthLevelDiagnosticsP,1) % loop over HOT and COLD definitions %!AG!
        [~,depthLevelDiagnosticsIdx(idP,idD)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelDiagnosticsP(idP,idD)));  %index of closest value to mean LAB depth further up in the plate
    end
end
UMdepthIdx          = depthLevelDiagnosticsIdx(:,1);
MMdepthIdx          = depthLevelDiagnosticsIdx(:,2);
LMdepthIdx          = depthLevelDiagnosticsIdx(:,3);
UMdepthP            = depthLevelDiagnosticsP(:,1); %!AG!  
MMdepthP            = depthLevelDiagnosticsP(:,2); %!AG!  
LMdepthP            = depthLevelDiagnosticsP(:,3); %!AG! 


%extracting data
for idD=1:size(depthLevelDiagnosticsP,2) % loop over UM, MM, LM
    for idP=1:size(depthLevelDiagnosticsP,1) % loop over HOT and COLD definitions %!AG!    
        %horizontal rms values
        if strcmp(GRID.Type,'yinyang')
            VHmax(idP,idD)   	= max(max(abs( [VH_3D(:,:,depthLevelDiagnosticsIdx(idP,idD)),VH_3Dyang(:,:,depthLevelDiagnosticsIdx(idP,idD))] )));
            VHmean(idP,idD)   	= mean2(abs( [VH_3D(:,:,depthLevelDiagnosticsIdx(idP,idD)),VH_3Dyang(:,:,depthLevelDiagnosticsIdx(idP,idD))] ));
        else
            VHmax(idP,idD)   	= max(max(abs( VH_3D(:,:,depthLevelDiagnosticsIdx(idP,idD)) )));
            VHmean(idP,idD)   	= mean2(abs( VH_3D(:,:,depthLevelDiagnosticsIdx(idP,idD)) ));
        end
    end
end

%writing it to output variables
MANTLE.VHmaxUM                 	= VHmax(:,1);    	%[plotting dim]
MANTLE.VHmaxMM                 	= VHmax(:,2);    	%[plotting dim]
MANTLE.VHmaxLM                 	= VHmax(:,3);    	%[plotting dim]
MANTLE.VHmeanUM               	= VHmean(:,1);    	%[plotting dim]
MANTLE.VHmeanMM                	= VHmean(:,2);    	%[plotting dim]
MANTLE.VHmeanLM               	= VHmean(:,3);    	%[plotting dim]
%!AG! 

%% CALCULATE RESIDUAL TEMPERATURE FIELD
%initiation
Tres3D = zeros(size(T_3D));

if MD.TresidualMode==1 %Standard horizontal residual
    SWITCH.customFieldTask  = 'Horizontal residual';
elseif MD.TresidualMode==2 %Horizontal-band residual
    SWITCH.customFieldTask  = 'Horizontal-band residual';
elseif MD.TresidualMode==3 %Global residual
    SWITCH.customFieldTask  = 'Global residual';
elseif MD.TresidualMode==4 %Regional residual
    SWITCH.customFieldTask  = 'Regional residual';
end
[VAR] = f_makeField(1,FILE,GRID,SETUP,SWITCH,PLOT,1,1);
%convert variables
if strcmp(GRID.Type,'yinyang')
    Tres3D              = VAR.var2d;
    Tres3Dyang          = VAR.var2d_yang;
    %ADD horizMean etc here for yin and yang ..................
    warning('from here on not implemented for yinyang grid...')
    return
    
    horizMeanTres       = VAR.meanTres; %might be regional values
    horizMinTres        = VAR.minTres;
    horizMaxTres        = VAR.maxTres;
else
    if strcmp(GRID.Dim,'2-D')
        Tres3D(:,1,:) 	= VAR.var2d;
    else %strcmp(GRID.Dim,'3-D')
        Tres3D          = VAR.var2d;
    end
    horizMeanTres       = VAR.meanTres; %might be regional values
    horizMinTres        = VAR.minTres;
    horizMaxTres        = VAR.maxTres;
end
clearvars VAR

%% CALCULATE RESIDUAL VZ FIELD 
%!AG!
%initiation
VZres3D = zeros(size(T_3D));
SWITCH.customFieldTask  = 'VZ Horizontal residual'; % additional field in f__makeField
[VAR] = f_makeField(1,FILE,GRID,SETUP,SWITCH,PLOT,1,1);
%convert variables
if strcmp(GRID.Type,'yinyang')
    VZres3D             = VAR.var2d;
    VZres3Dyang         = VAR.var2d_yang;
    %ADD horizMean etc here for yin and yang ..................
    warning('from here on not implemented for yinyang grid...')
    return
    
    horizMeanVZres    	= VAR.meanVZres; %might be regional values
    horizMinVZres   	= VAR.minVZres;
    horizMaxVZres    	= VAR.maxVZres;
else
    if strcmp(GRID.Dim,'2-D')
        VZres3D(:,1,:) 	= VAR.var2d;
    else %strcmp(GRID.Dim,'3-D')
        VZres3D        	= VAR.var2d;
    end
    horizMeanVZres   	= VAR.meanVZres; %might be regional values
    horizMinVZres       = VAR.minVZres;
    horizMaxVZres    	= VAR.maxVZres;
end
clearvars VAR

%% CALCULATE RESIDUAL RADIAL HEAT ADVECTION FIELD
%!AG! 
% initiation
VZresTres3D = zeros(size(Tres3D));
horizMeanVZresTres  = VZresTres3D;
horizMinVZresTres   = VZresTres3D;
horizMaxVZresTres   = VZresTres3D;

% calculation
% vzTres = vzres.*abs(Tres)
%                 abs(Tres): to exclude downwellings (both VZres and Tres -ive)
VZresTres3D         = VZres3D.*abs(Tres3D);
horizMeanVZresTres  = horizMeanVZres.*abs(horizMeanTres);
horizMinVZresTres   = horizMinVZres .*abs(horizMinTres);
horizMaxVZresTres   = horizMeanVZres.*abs(horizMaxTres);

sec_in_yr      = 31556926; % seconds in a year !AG! to do -- different!
if logical(0)
    figure(33); clf
    x(:,:) = GRID.X_3Dp(:,1,:);
    z(:,:) = GRID.Z_3Dp(:,1,:);
    p(:,:) = VZresTres3D(:,1,:);
    pcolor(x,z,p*sec_in_yr)
    colorbar
    axis ij
end

%% FIND UP- AND DOWNWELLINGS
thrVZ                           = 1/100 *VZmax;       	%define threshold for up/down welling <<<<<<<<<<<<<<<<<<
thrVZup                         = 1/100 *VZupmax;       %define threshold for upwelling <<<<<<<<<<<<<<<<<<
thrVZdown                       = 1/100 *VZdownmax;   	%define threshold for downwelling <<<<<<<<<<<<<<<<<<
upWelling                       = zeros(size(VZ_3D));  	%nothing
downWelling                     = zeros(size(VZ_3D));   %nothing
upWellingAbsolute             	= zeros(size(VZ_3D));  	%nothing
downWellingAbsolute            	= zeros(size(VZ_3D));   %nothing
upWelling(VZ_3D>thrVZup)        = 1;                  	%upwelling (exceeding threshold)
downWelling(VZ_3D<-thrVZdown)	= 1;                    %downwelling (exceeding threshold)
upWellingAbsolute(VZ_3D>0)      = 1;                  	%upwelling (absolute)
downWellingAbsolute(VZ_3D<0)    = 1;                  	%downwelling (absolute)
if strcmp(GRID.Type,'yinyang')
    upWelling_yang                          = zeros(size(VZ_3Dyang)); 	%nothing
    downWelling_yang                        = zeros(size(VZ_3Dyang)); 	%nothing
    upWellingAbsolute_yang              	= zeros(size(VZ_3Dyang)); 	%nothing
    downWellingAbsolute_yang             	= zeros(size(VZ_3Dyang)); 	%nothing
    upWelling_yang(VZ_3Dyang>thrVZup)       = 1;                        %upwelling (exceeding threshold)
    downWelling_yang(VZ_3Dyang<-thrVZdown)  = 1;                        %downwelling (exceeding threshold)
    upWellingAbsolute_yang(VZ_3Dyang>0)   	= 1;                        %upwelling (absolute)
    downWellingAbsolute_yang(VZ_3Dyang<0) 	= 1;                        %downwelling (absolute)
end

%% FIND ACTIVE PLUMES 
% (1) thermally, like in Labrosse (EPSL,2002); 
% (2) thermal and thermal*dynamical, using percentiles (and absolute minimum thresholds !AG!
%Initialising
Plumes3D = zeros(size(Tres3D));
thrHot  = Plumes3D;
thrCold = Plumes3D;
%make thresholds 
sec_in_yr      = 31556926; % seconds in a year
thrHot_T       = 0;
thrHot_VZT     = 0;
thrCold_T      = 0;
thrCold_VZT    = 0;
thrHot_VZTp    = 90;       % make into PLOT.VZTperHOT
thrHot_Tp      = 85;       % make into PLOT.TperHOT
thrCold_VZTp   = 85;       % make into PLOT.VZTperCOLD
thrCold_Tp     = 90;       % make into PLOT.TperCOLD
thr_Tres_abs   = 50;                                % CHECK 
thr_VZTres_abs = 0.01/sec_in_yr * thr_Tres_abs;     % CHECK 
w3D = GRID.cellVolume/(sum(GRID.cellVolume(:))*nb); 

if strcmp(GRID.Type,'yinyang'); Plumes3Dyang = Plumes3D; end

% HOT PLUMES
if MD.activeUpwellingDefinition==1
    %start hot plume diagnostics
    thrHot                    	= horizMeanTres + PLOT.pHot*(horizMaxTres-horizMeanTres);	%after Labrosse (EPSL,2002)
    
    Plumes3D(Tres3D>thrHot)   	= 1;                    %hot plumes
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang(Tres3Dyang>thrHot_yang)    = 1;  	%hot plumes 
    end
elseif MD.activeUpwellingDefinition==2 %!AG!
    % using weighted percentiles (on cell volume)
    thrHot_T   = f_wprctile(reshape(Tres3D,[],1),thrHot_Tp,reshape(w3D,[],1));
    thrHot_VZT = f_wprctile(reshape(VZresTres3D,[],1),thrHot_VZTp,reshape(w3D,[],1));
    
    Plumes3D(Tres3D>thrHot_T  & Tres3D>thr_Tres_abs & ...
             VZresTres3D>thrHot_VZT & VZresTres3D>thr_VZTres_abs) = 1;
end

% COLD PLUMES
if MD.activeDownwellingDefinition==1
    %start cold plume diagnostics
    thrCold                 	= horizMeanTres + PLOT.pCold.*(horizMinTres-horizMeanTres);	%after Labrosse (EPSL,2002)
    
    Plumes3D(Tres3D<thrCold & Tres3D<(-thr_Tres_abs))  	= -1;                   %cold plumes
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang(Tres3Dyang<thrCold_yang)  	= -1;  	%cold plumes
    end
elseif MD.activeDownwellingDefinition==2 %!AG!
    % using weighted percentiles (on cell volume)
    thrCold_T   = f_wprctile(reshape((-Tres3D),[],1),thrCold_Tp,reshape(w3D,[],1));        % use -ive because downwellings
    thrCold_VZT = f_wprctile(reshape((-VZresTres3D),[],1),thrCold_VZTp,reshape(w3D,[],1)); % use -ive because downwellings
    
    Plumes3D((-Tres3D)>thrCold_T  & (-Tres3D)>thr_Tres_abs & ...                 % use -ive because downwellings
             (-VZresTres3D)>thrCold_VZT & (-VZresTres3D)>thr_VZTres_abs) = -1;    % use -ive because downwellings
end

clearvars dummy dummyX dummyY
clearvars horizMeanT horizMeanTres horizMax horizMin thrHot thrCold
clearvars thrHot_T thrHot_VZT thrCold_T thrCold_VZT thrHot_Tp thrHot_VZTp thrCold_Tp thrCold_VZTp ... 
          thr_Tres_abs thr_VZTres_abs w3D %!AG!

if logical(0)
    figure(3)
    hold on
    x(:,:) = GRID.X_3Dp(:,1,:);
    z(:,:) = GRID.Z_3Dp(:,1,:);
    p(:,:) = Plumes3D(:,1,:);
    phot = zeros(size(p));
    phot(p==1) = 1;
    phot = bwmorph(phot,'thicken',1);
    phot = bwmorph(phot,'bridge');
    
    pcold = zeros(size(p));
    pcold(p==-1) = 1;
    %    pcold = bwmorph(pcold,'thicken',2);
    pcold = bwmorph(pcold,'bridge');
    
    p2 = phot-pcold;
    contourf(x,z,p2)
    colorbar
    axis ij
    figure(1)
end

%% FIND Continents
%critical variables
thrCont                     = 0.5;  %composional threshold for continental material (0<1) <<<<<<<<<<<<<<<<<<<<<
%initialising
CONT3D = zeros(size(CC_3D));
if strcmp(GRID.Type,'yinyang'); CONT3Dyang = CONT3D; end
%start diagnostics
CONT3D(CC_3D>thrCont)   	= 1;            	%continental material

%% FIND LLSVPs
%critical variables
thrPrm                      = 0.5;  %composional threshold for primordial material (0<1) <<<<<<<<<<<<<<<<<<<<<
%initialising
LLSVP3D = zeros(size(PRM_3D));
if strcmp(GRID.Type,'yinyang'); LLSVP3Dyang = LLSVP3D; end
%start diagnostics
%LLSVP3D(PRM_3D>thrPrm)   	= 1;            	%llsvp material
%!AG!
thrBs                       = 0.5;
thrT                        = 3333;
LLSVP3D(BS_3D>thrBs & T_3D>thrT)          = 1;                    %llsvp material
%!AG! re-define LLSVP's??
%C_bas > 0.75
%T_cell > (3000 K + T_CMB)/2 
% ^ Schierjott et al 2020

%% ACCOUNT FOR PERIODIC BOUNDARIES
%check for plumes crossing periodic boundaries:
add_factor          = MD.addFactor;      %percent of total x-extent that is added at x=1 and x=end
%add some slices at x=1 and x=end
PlumesOrig          = Plumes3D;
extentXgrid         = false;
extentYgrid         = false;
nx_orig             = size(PlumesOrig,1);
ny_orig             = size(PlumesOrig,2);
nx_add              = nx_orig*add_factor;
ny_add              = ny_orig*add_factor;
if strcmp(GRID.Type,'yinyang')
    extentYgrid     = true;
elseif strcmp(GRID.Type,'Cartesian')
    if strcmp(GRID.Dim,'2-D')
        if size(Plumes3D,1)>1 && size(Plumes3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
        if size(Plumes3D,1)==1 && size(Plumes3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
    elseif strcmp(GRID.Dim,'3-D')
        if strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
        if strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
        warning(['MANTLE DIAGNOSTICS: Grid type ',GRID.Type,' in 3-D not tested yet!'])
    end
elseif strcmp(GRID.Type,'spherical2D')
    if size(Plumes3D,1)>1 && size(Plumes3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
    if size(Plumes3D,1)==1 && size(Plumes3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
else
    error(['Grid type ',GRID.Type,' not found!'])
end
if extentXgrid
    X_3DpExt        = [GRID.X_3Dp((end-nx_add+1):end,:,:); GRID.X_3Dp; GRID.X_3Dp(1:nx_add,:,:)];
    Y_3DpExt        = [GRID.Y_3Dp((end-nx_add+1):end,:,:); GRID.Y_3Dp; GRID.Y_3Dp(1:nx_add,:,:)];
    Z_3DpExt        = [GRID.Z_3Dp((end-nx_add+1):end,:,:); GRID.Z_3Dp; GRID.Z_3Dp(1:nx_add,:,:)];
    Plumes3D        = [Plumes3D((end-nx_add+1):end,:,:); Plumes3D; Plumes3D(1:nx_add,:,:)];
    CONT3D          = [CONT3D((end-nx_add+1):end,:,:); CONT3D; CONT3D(1:nx_add,:,:)];
    LLSVP3D         = [LLSVP3D((end-nx_add+1):end,:,:); LLSVP3D; LLSVP3D(1:nx_add,:,:)];
    upWelling       = [upWelling((end-nx_add+1):end,:,:); upWelling; upWelling(1:nx_add,:,:)];
    downWelling     = [downWelling((end-nx_add+1):end,:,:); downWelling; downWelling(1:nx_add,:,:)];
    upWellingAbsolute       = [upWellingAbsolute((end-nx_add+1):end,:,:); upWellingAbsolute; upWellingAbsolute(1:nx_add,:,:)];
    downWellingAbsolute     = [downWellingAbsolute((end-nx_add+1):end,:,:); downWellingAbsolute; downWellingAbsolute(1:nx_add,:,:)];
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang 	= [Plumes3Dyang((end-nx_add+1):end,:,:); Plumes3Dyang; Plumes3Dyang(1:nx_add,:,:)];
        LLSVP3Dyang 	= [LLSVP3Dyang((end-nx_add+1):end,:,:); LLSVP3Dyang; LLSVP3Dyang(1:nx_add,:,:)];
        upWelling_yang 	= [upWelling_yang((end-nx_add+1):end,:,:); upWelling_yang; upWelling_yang(1:nx_add,:,:)];
        downWelling_yang= [downWelling_yang((end-nx_add+1):end,:,:); downWelling_yang; downWelling_yang(1:nx_add,:,:)];
        upWellingAbsolute_yang 	= [upWellingAbsolute_yang((end-nx_add+1):end,:,:); upWellingAbsolute_yang; upWellingAbsolute_yang(1:nx_add,:,:)];
        downWellingAbsolute_yang= [downWellingAbsolute_yang((end-nx_add+1):end,:,:); downWellingAbsolute_yang; downWellingAbsolute_yang(1:nx_add,:,:)];
    end
else
    X_3DpExt        = GRID.X_3Dp;
    Y_3DpExt        = GRID.X_3Dp;
    Z_3DpExt        = GRID.X_3Dp;
end
if extentYgrid
    X_3DpExt     	= [X_3DpExt(:,(end-ny_add+1):end,:), X_3DpExt, X_3DpExt(:,1:ny_add,:)];
    Y_3DpExt     	= [Y_3DpExt(:,(end-ny_add+1):end,:), Y_3DpExt, Y_3DpExt(:,1:ny_add,:)];
    Z_3DpExt     	= [Z_3DpExt(:,(end-ny_add+1):end,:), Z_3DpExt, Z_3DpExt(:,1:ny_add,:)];
    Plumes3D     	= [Plumes3D(:,(end-ny_add+1):end,:), Plumes3D, Plumes3D(:,1:ny_add,:)];
    CONT3D          = [CONT3D(:,(end-ny_add+1):end,:), CONT3D, CONT3D(:,1:ny_add,:)];
    LLSVP3D     	= [LLSVP3D(:,(end-ny_add+1):end,:), LLSVP3D, LLSVP3D(:,1:ny_add,:)];
    upWelling       = [upWelling(:,(end-ny_add+1):end,:), upWelling, upWelling(:,1:ny_add,:)];
    downWelling     = [downWelling(:,(end-ny_add+1):end,:), downWelling, downWelling(:,1:ny_add,:)];
    upWellingAbsolute       = [upWellingAbsolute(:,(end-ny_add+1):end,:), upWellingAbsolute, upWellingAbsolute(:,1:ny_add,:)];
    downWellingAbsolute     = [downWellingAbsolute(:,(end-ny_add+1):end,:), downWellingAbsolute, downWellingAbsolute(:,1:ny_add,:)];
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang 	= [Plumes3Dyang(:,(end-ny_add+1):end,:), Plumes3Dyang, Plumes3Dyang(:,1:ny_add,:)];
        LLSVP3Dyang 	= [LLSVP3Dyang(:,(end-ny_add+1):end,:), LLSVP3Dyang, LLSVP3Dyang(:,1:ny_add,:)];
        upWelling_yang	= [upWelling_yang(:,(end-ny_add+1):end,:), upWelling_yang, upWelling_yang(:,1:ny_add,:)];
        downWelling_yang= [downWelling_yang(:,(end-ny_add+1):end,:), downWelling_yang, downWelling_yang(:,1:ny_add,:)];
        upWellingAbsolute_yang	= [upWellingAbsolute_yang(:,(end-ny_add+1):end,:), upWellingAbsolute_yang, upWellingAbsolute_yang(:,1:ny_add,:)];
        downWellingAbsolute_yang= [downWellingAbsolute_yang(:,(end-ny_add+1):end,:), downWellingAbsolute_yang, downWellingAbsolute_yang(:,1:ny_add,:)];
    end
end

for ib=1:nb
    Plumes3Dh = zeros(size(Plumes3D)); Plumes3Dc = Plumes3Dh;
    CONT3Dx = zeros(size(Plumes3D));
    LLSVP3Dx = zeros(size(Plumes3D));
    if ib==1
        Plumes3Dh(Plumes3D==1)  = 1;
        Plumes3Dc(Plumes3D==-1) = 1;
        CONT3Dx                 = CONT3D;
        LLSVP3Dx                = LLSVP3D;
    elseif ib==2
        Plumes3Dh(Plumes3Dyang==1)  = 1;
        Plumes3Dc(Plumes3Dyang==-1) = 1;
        CONT3Dx                 = CONT3Dyang;
        LLSVP3Dx                = LLSVP3Dyang;
    end
    %% OPTIMISE DATA I
    if logical(1)
        if strcmp(GRID.Dim,'2-D')
            %Plumes
            dummyh(:,:)         = Plumes3Dh(:,1,:); %make 2-dimensional
            dummyc(:,:)         = Plumes3Dc(:,1,:); %make 2-dimensional
            %thicken
%             dummyh              = bwmorph(dummyh,'thicken',1);
%             dummyc              = bwmorph(dummyc,'thicken',1);
            %bridge
            dummyh              = bwmorph(dummyh,'bridge');
            dummyc              = bwmorph(dummyc,'bridge');
            Plumes3Dh(:,1,:)    = dummyh(:,:);
            Plumes3Dc(:,1,:)    = dummyc(:,:);
            clearvars dummyh dummyc
            
            for iF=1:2
                if iF==1 %Continents
                    dummy(:,:)          = CONT3Dx(:,1,:); %make 2-dimensional
                elseif iF==2 %LLSVPs
                    dummy(:,:)          = LLSVP3Dx(:,1,:); %make 2-dimensional
                else
                    error('check here')
                end
                %thicken
                %             dummy              = bwmorph(dummy,'thicken',1);
                %bridge
                dummy               = bwmorph(dummy,'bridge');
                %fill isolated pixels
                dummy               = bwmorph(dummy,'fill');
                %fill holes
                dummy               = imfill(dummy,'holes');   %<<<<<<<<<if it is filling large holes, the volume diagnostics will be inaccurate!
                if iF==1 %Continents
                    CONT3Dx(:,1,:)    	= dummy(:,:);
                elseif iF==2 %LLSVPs
                    LLSVP3Dx(:,1,:)    	= dummy(:,:);
                end
                clearvars dummy
            end
            
        elseif strcmp(GRID.Dim,'3-D')
            if SWITCH.Verbose; warning('Optimising plume data in 3-D not possible yet'); end
        end
    end
    
    if logical(0)
        figure(3)
        x(:,:) = X_3DpExt(:,1,:);
        z(:,:) = Z_3DpExt(:,1,:);
        p(:,:) = Plumes3Dh(:,1,:);
        contourf(x,z,p)
        colorbar
        axis ij
        figure(1)
    end
   
    %% CHECK CONNECTIVITY FOR AREA SIZES
    CCph            = bwconncomp(Plumes3Dh);     	%HOT PLUMES
    numPixelsPh     = cellfun(@numel,CCph.PixelIdxList);    %number of connected pixels in each connected area
    CCpc            = bwconncomp(Plumes3Dc);     	%COLD PLUMES
    numPixelsPc     = cellfun(@numel,CCpc.PixelIdxList);    %number of connected pixels in each connected area
    CCuw            = bwconncomp(upWelling);        %UPWELLINGS
    numPixelsUW     = cellfun(@numel,CCuw.PixelIdxList);    %number of connected pixels in each connected area
    CCdw            = bwconncomp(downWelling);      %DOWNWELLINGS
    numPixelsDW     = cellfun(@numel,CCdw.PixelIdxList);    %number of connected pixels in each connected area
    CCcont          = bwconncomp(CONT3Dx);       	%Continents
    numPixelsCont	= cellfun(@numel,CCcont.PixelIdxList); %number of connected pixels in each connected area
    CCllsvp         = bwconncomp(LLSVP3Dx);       	%LLSVPs
    numPixelsLLSVP	= cellfun(@numel,CCllsvp.PixelIdxList); %number of connected pixels in each connected area
    
    %% REMOVING SMALL ANOMALIES
    removeSmallAnomalies    = logical(1);
    plumesBhot  = Plumes3Dh;  	%hot plumes
    plumesBcold = Plumes3Dc;  	%cold plumes
    if removeSmallAnomalies
        %Plumes
        % *********************** small area criterion <<<<<<<<<<<<<<<<
        % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
        sizeThreshold           = 1/5* size(Plumes3Dh,3);   %find areas larger than (1/5*nz) pixels
        if strcmp(GRID.Dim,'3-D')
            sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
        end
        indLargePh          	= find(numPixelsPh<sizeThreshold);
        indLargePc          	= find(numPixelsPc<sizeThreshold);
        % ***********************
        for i_ind=1:size(indLargePh,2)
            plumesBhot(CCph.PixelIdxList{indLargePh(i_ind)})    = 0;
        end
        for i_ind=1:size(indLargePc,2)
            plumesBcold(CCpc.PixelIdxList{indLargePc(i_ind)})   = 0;
        end
        
        for iF=1:2
            if iF==1 %Continents
                % *********************** small area criterion <<<<<<<<<<<<<<<<
                % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
                sizeThreshold           = 10* size(CONT3Dx,3);   %find areas larger than (1/5*nz) pixels
                if strcmp(GRID.Dim,'3-D')
                    sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
                end
                indLargeCont          	= find(numPixelsCont<sizeThreshold);
                % ***********************
                for i_ind=1:size(indLargeCont,2)
                    CONT3Dx(CCcont.PixelIdxList{indLargeCont(i_ind)}) = 0;
                end
            elseif iF==2 %LLSVPs
                % *********************** small area criterion <<<<<<<<<<<<<<<<
                % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
                sizeThreshold           = 2* size(LLSVP3Dx,3);   %find areas larger than (1/15*nz) pixels
                if strcmp(GRID.Dim,'3-D')
                    sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
                end
                indLargeLlsvp          	= find(numPixelsLLSVP<sizeThreshold);
                % ***********************
                for i_ind=1:size(indLargeLlsvp,2)
                    LLSVP3Dx(CCllsvp.PixelIdxList{indLargeLlsvp(i_ind)}) = 0;
                end
            end
        end
    end
    
    %% SETTING UP HOT AND COLD ANOMALIES after criterion
    if ib==1
        plumesHot             	= plumesBhot;   %hot plumes
        plumesCold              = plumesBcold;  %cold plumes
        CONT3D                  = CONT3Dx;      %continents
        LLSVP3D                 = LLSVP3Dx;     %llsvp
    elseif ib==2
        plumesHot_yang      	= plumesBhot;   %hot plumes
        plumesCold_yang     	= plumesBcold;	%cold plumes
        CONT3Dyang              = CONT3Dx;      %continents
        LLSVP3Dyang           	= LLSVP3Dx;     %llsvp
    end
    clearvars plumesBhot plumesBcold
end

%% DISCRIMINATE BETWEEN ACTIVE UPWELLINGS AND LLSVP'S !AG!
Plumes3D(ones(size(Plumes3D)) & LLSVP3D) = 0;
plumesHot(ones(size(plumesHot)) & LLSVP3D) = 0;

if logical(0)
    figure(3)
    x(:,:) = X_3DpExt(:,1,:);
    z(:,:) = Z_3DpExt(:,1,:);
%    p(:,:) = plumesHot(:,1,:);
     p(:,:) = LLSVP3D(:,1,:);
    contourf(x,z,p)
    colorbar
    axis ij
    figure(1)
end

%% CONNECTIVITY CHECK WITH CORRESPONDING BOUNDARY LAYER
% ...AND CHECK FOR PLUME EXTENSION THROUGHOUT UPPER-AND LOWER-DEPTH THRESHOLDS
% This uses PLOT.Z_3Dp, which is always flipped and actual depthCCPhot   = bwconncomp(plumesHot);
CCPhot      = bwconncomp(plumesHot);
CCPcold     = bwconncomp(plumesCold);
if strcmp(GRID.Type,'yinyang')
    CCPhot_yang   = bwconncomp(plumesHot_yang);
    CCPcold_yang  = bwconncomp(plumesCold_yang);
end
Ztop    = 0;
Zbot    = max(Z_3DpExt(:));
Zdiff   = Zbot-Ztop;
if Zdiff<=0; error('Zdiff should be >0'); end
for ib=1:nb
    for plume_kind=1:2
        if plume_kind==1  %HOT
            if ib==1
                CC_kind = CCPhot; 
            else
                CC_kind = CCPhot_yang; 
            end
            ZtopThreshold  	= Ztop + Zdiff*MD.upperDepthThresholdHot;
            ZbotThreshold  	= Ztop + Zdiff*MD.lowerDepthThresholdHot;
            dzThreshold     = Zdiff*MD.dzThresholdHot; %!AG!
        elseif plume_kind==2  %COLD
            if ib==1
                CC_kind = CCPcold;
            elseif ib==2
                CC_kind = CCPcold_yang;
            end
            ZtopThreshold  	= Ztop + Zdiff*MD.upperDepthThresholdCold;
            ZbotThreshold 	= Ztop + Zdiff*MD.lowerDepthThresholdCold;
            dzThreshold     = Zdiff*MD.dzThresholdCold; %!AG!
        end
        
        if (MD.DepthThresholdMode == 1)
            for i_area=1:size(CC_kind.PixelIdxList,2)  %loop all connected areas
                surf_connection = false;
                bot_connection = false;
                for ip=1:size(CC_kind.PixelIdxList{1,i_area},1)  %loop all pixels of an area
                    i_pixel = CC_kind.PixelIdxList{1,i_area}(ip,1);
                    %[x,y,z] = ind2sub(size(PLUMES),i_pixel);
                    if Z_3DpExt(i_pixel)<=ZtopThreshold %is connected to the top (0)
                        surf_connection = true;
                    end
                    if Z_3DpExt(i_pixel)>=ZbotThreshold %is connected to the bottom (1)
                        bot_connection = true;
                    end
                    if surf_connection && bot_connection
                        break   %this is a plume exceeding top to bottom levels! check next.
                    end
                end
                
                if surf_connection && bot_connection %of defined depth thresholds
                %keep area
                else %remove area: this is no plume
                    if plume_kind==1  %HOT
                        if ib==1
                            plumesHot(CC_kind.PixelIdxList{i_area})         = 0;
                        elseif ib==2
                            plumesHot_yang(CC_kind.PixelIdxList{i_area})    = 0;
                        end
                    elseif plume_kind==2  %COLD
                        if ib==1
                            plumesCold(CC_kind.PixelIdxList{i_area})        = 0;
                        elseif ib==2
                            plumesCold_yang(CC_kind.PixelIdxList{i_area})   = 0;
                        end
                    end
                end
            end
        elseif(MD.DepthThresholdMode == 2) %!AG!  
            % Find top; bottom depths of each connected area, 
            % if (ztop-zbot)<dz_threshold: remove area 
            for i_area=1:size(CC_kind.PixelIdxList,2)  %loop all connected areas
                dz_connection = false;
                zMin_area = 0;
                zMax_area = 0; 
            
                i_pixels = CC_kind.PixelIdxList{1,i_area}(:,1); % i_pixel array
                %[x,y,z] = ind2sub(size(PLUMES),i_pixel);
                zMin_area = min(Z_3DpExt(i_pixels(:))); % minimal Z value in area
                zMax_area = max(Z_3DpExt(i_pixels(:))); % maximal Z value in area
                if((zMax_area-zMin_area)>dzThreshold)
                    dz_connection = true;
                end
                
                %!AG! add part where slabs in the UM are still detected,
                %even though they do not pass the dz_threshold test
                if (plume_kind == 2)
                    if (zMin_area<= ZtopThreshold && zMax_area >= ZbotThreshold)
                        dz_connection = true;
                    end
                end
            
                if dz_connection 
                    %keep area
                else
                    %remove area: this is no plume
                    if plume_kind==1  %HOT
                        if ib==1
                            plumesHot(CC_kind.PixelIdxList{i_area})         = 0;
                        elseif ib==2
                            plumesHot_yang(CC_kind.PixelIdxList{i_area})    = 0;
                        end
                    elseif plume_kind==2  %COLD
                        if ib==1
                            plumesCold(CC_kind.PixelIdxList{i_area})        = 0;
                        elseif ib==2
                            plumesCold_yang(CC_kind.PixelIdxList{i_area})   = 0;
                        end
                    end
                end
            end
        end
    end
end


%% GET CENTROIDS FOR EXTRACTED ANOMALIES
if ~strcmp(GRID.Type,'yinyang') && ~strcmp(GRID.Dim,'3-D')
    for iF=1:2
        if iF==1 %continents
            CCanom   	= bwconncomp(CONT3D);
            % CCllsvp2_yang   	= bwconncomp(CONT3Dyang);
        elseif iF==2 && ~(sum(sum(LLSVP3D))==0) %llsvp
            CCanom   	= bwconncomp(LLSVP3D) %!AG! somehow this is non-zero even when there are no LLSVP detected!
                                              %!AG! that's why the added condition ~(sum(sum(LLSVP3D))==0)
            % CCllsvp2_yang   	= bwconncomp(LLSVP3Dyang);
        end
         CentroidAnomIdx = zeros(size(CCanom.PixelIdxList,2),2);
        for i_area=1:size(CCanom.PixelIdxList,2)  %loop all connected areas
            areaCurrent(:,:)    	= false(size(LLSVP3D));
            areaCurrent(CCanom.PixelIdxList{i_area})    = true;
            %     areaCurrent     = bwareafilt(areaCurrent, 1);
            
            props                   = regionprops(areaCurrent, 'Centroid');
            dummy(i_area,:)         = round(props.Centroid); %find closest index
            CentroidAnomIdx(i_area,1)	= min(max(1,dummy(i_area,2)),size(LLSVP3D,1)); %need to flip x and z values, and limit to max/min indices
            CentroidAnomIdx(i_area,2)	= min(max(1,dummy(i_area,1)),size(LLSVP3D,3));
            
            % binaryImage = bwareafilt(binaryImage, 1); % Extract largest blob ONLY.
            % props = regionprops(binaryImage, 'Centroid');
            % centroid = props.Centroid;
        end
        if iF==1 && ~(sum(sum(CONT3D))==0) %continents
            %convert to actual plot position
            dummy2(1,:)             = Z_3DpExt(1,1,CentroidAnomIdx(:,2));

            CentroidCONT       = [X_3DpExt(CentroidAnomIdx(:,1),1,1), dummy2'];
            %remove dublicated entries (due to periodic side boundaries)
            CentroidCONT       = unique(CentroidCONT,'rows');
        elseif iF==2 && ~(sum(sum(LLSVP3D))==0) %llsvps
            %convert to actual plot position
            dummy2(1,:)             = Z_3DpExt(1,1,CentroidAnomIdx(:,2));

            CentroidLLSVP       = [X_3DpExt(CentroidAnomIdx(:,1),1,1), dummy2'];
            %remove dublicated entries (due to periodic side boundaries)
            CentroidLLSVP       = unique(CentroidLLSVP,'rows');
        end
        clearvars dummy dummy2 areaCurrent

        if logical(0)
            figure(3)
            x(:,:) = X_3DpExt(:,1,:);
            z(:,:) = Z_3DpExt(:,1,:);
            if iF==1 %continents
                p(:,:) = CONT3D(:,1,:);
                centroidDummy = CentroidCONT;
            else
                p(:,:) = LLSVP3D(:,1,:);
                centroidDummy = CentroidLLSVP;
            end
            contourf(x,z,p)
            hold on
            for iiii=1:size(CentroidAnomIdx,1)
                plot(X_3DpExt(CentroidAnomIdx(iiii,1),1,1),Z_3DpExt(1,1,CentroidAnomIdx(iiii,2)),'r+')
            end
            plot(centroidDummy(:,1),centroidDummy(:,2),'gx')
            colorbar
            axis ij
            figure(1)
        end
    end
else
    
    %add for yinyang!
    
end

%% REMOVING GRID EXTENSION (BACK TO ORIGINAL SIZE)
if extentXgrid
    Plumes3D            = Plumes3D((nx_add+1):end-nx_add,:,:);
    plumesHot           = plumesHot((nx_add+1):end-nx_add,:,:);
    plumesCold          = plumesCold((nx_add+1):end-nx_add,:,:);
    upWelling           = upWelling((nx_add+1):end-nx_add,:,:);
    downWelling         = downWelling((nx_add+1):end-nx_add,:,:);
    upWellingAbsolute  	= upWellingAbsolute((nx_add+1):end-nx_add,:,:);
    downWellingAbsolute	= downWellingAbsolute((nx_add+1):end-nx_add,:,:);
    CONT3D              = CONT3D((nx_add+1):end-nx_add,:,:);
    LLSVP3D             = LLSVP3D((nx_add+1):end-nx_add,:,:);
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang           = Plumes3Dyang((nx_add+1):end-nx_add,:,:);
        plumesHot_yang          = plumesHot_yang((nx_add+1):end-nx_add,:,:);
        plumesCold_yang         = plumesCold_yang((nx_add+1):end-nx_add,:,:);
        upWelling_yang          = upWelling_yang((nx_add+1):end-nx_add,:,:);
        downWelling_yang        = downWelling_yang((nx_add+1):end-nx_add,:,:);
        upWellingAbsolute_yang	= upWellingAbsolute_yang((nx_add+1):end-nx_add,:,:);
        downWellingAbsolute_yang= downWellingAbsolute_yang((nx_add+1):end-nx_add,:,:);
        CONT3Dyang             = CONT3Dyang((nx_add+1):end-nx_add,:,:);
        LLSVP3Dyang            = LLSVP3Dyang((nx_add+1):end-nx_add,:,:);
    end
end
if extentYgrid
    Plumes3D            = Plumes3D(:,(ny_add+1):end-ny_add,:);
    plumesHot           = plumesHot(:,(ny_add+1):end-ny_add,:);
    plumesCold          = plumesCold(:,(ny_add+1):end-ny_add,:);
    upWelling           = upWelling(:,(ny_add+1):end-ny_add,:);
    downWelling         = downWelling(:,(ny_add+1):end-ny_add,:);
    upWellingAbsolute  	= upWellingAbsolute(:,(ny_add+1):end-ny_add,:);
    downWellingAbsolute	= downWellingAbsolute(:,(ny_add+1):end-ny_add,:);
    CONT3D              = CONT3D(:,(ny_add+1):end-ny_add,:);
    LLSVP3D             = LLSVP3D(:,(ny_add+1):end-ny_add,:);
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang           = Plumes3Dyang(:,(ny_add+1):end-ny_add,:);
        plumesHot_yang          = plumesHot_yang(:,(ny_add+1):end-ny_add,:);
        plumesCold_yang         = plumesCold_yang(:,(ny_add+1):end-ny_add,:);
        upWelling_yang          = upWelling_yang(:,(ny_add+1):end-ny_add,:);
        downWelling_yang        = downWelling_yang(:,(ny_add+1):end-ny_add,:);
        upWellingAbsolute_yang 	= upWellingAbsolute_yang(:,(ny_add+1):end-ny_add,:);
        downWellingAbsolute_yang= downWellingAbsolute_yang(:,(ny_add+1):end-ny_add,:);
        CONT3Dyang             = CONT3Dyang(:,(ny_add+1):end-ny_add,:);
        LLSVP3Dyang            = LLSVP3Dyang(:,(ny_add+1):end-ny_add,:);
    end
end
clearvars X_3DpExt Y_3DpExt Z_3DpExt


%% DESCRIMINATE BETWEEN ACTIVE AND PASSIVE UP-/DOWNWELLING
passiveUpwelling    = upWelling-plumesHot;      passiveUpwelling(passiveUpwelling<0) = 0;
passiveDownwelling 	= downWelling-plumesCold;   passiveDownwelling(passiveDownwelling<0) = 0;
%!AG!
activeUpwelling     = zeros(size(upWelling));   activeUpwelling(ones(size(upWelling)) & plumesHot) = 1;
activeDownwelling       = zeros(size(downWelling)); activeDownwelling(ones(size(downWelling)) & plumesCold) = 1;
%!AG!
%activeUpwelling     = zeros(size(upWelling));   activeUpwelling(upWellingAbsolute & plumesHot) = 1;
%activeDownwelling	= zeros(size(downWelling)); activeDownwelling(downWellingAbsolute & plumesCold) = 1;
if strcmp(GRID.Type,'yinyang')
    passiveUpwelling_yang 	= upWelling_yang-plumesHot_yang; 	passiveUpwelling_yang(passiveUpwelling_yang<0) = 0;
    passiveDownwelling_yang	= downWelling_yang-plumesCold_yang;	passiveDownwelling_yang(passiveDownwelling_yang<0) = 0;
    activeUpwelling_yang    = zeros(size(upWelling_yang));      activeUpwelling_yang(upWellingAbsolute_yang & plumesHot_yang) = 1;
    activeDownwelling_yang	= zeros(size(downWelling_yang));    activeDownwelling_yang(downWellingAbsolute_yang & plumesCold_yang) = 1;
end

%% LLSVP DEFINITION
MANTLE.llsvp = LLSVP3D


%% PLUME DEFINITION
if MD.plumeDefinition==1 %defined solely by temperature
%     plumesHot       = plumesHot;
%     plumesCold      = plumesCold;
elseif MD.plumeDefinition==2 %defined by temperature and vertical velocity
    plumesHot       = activeUpwelling;
    plumesCold      = activeDownwelling;
end


if ~strcmp(GRID.Type,'yinyang')
    %% DETAILS OF FOUND ANOMALIES I: number and lateral mobility at certain depths
    %derive horizontal slices
    depthLevelCONTslices    = 0.99;                         %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<<<<
    depthLevelLLSVPslices   = 0.01;                         %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<<<<
    depthLevelSlicesIdx     = [UMdepthIdx,MMdepthIdx,LMdepthIdx];
    depthLevelSlicesPL      = [UMdepthP,MMdepthP,LMdepthP];
    
    %PLUMES
    %derive horizontal plume slices
    plumesHotSlice      = plumesHot(:,:,depthLevelSlicesIdx(DepthIdxHot,:)); %!AG! aren't these the same ?? 
    plumesColdSlice   	= plumesCold(:,:,depthLevelSlicesIdx(DepthIdxCold,:));
    plumesHotSlice      = plumesHot(:,:,[UMdepthIdx(DepthIdxHot),MMdepthIdx(DepthIdxHot),LMdepthIdx(DepthIdxHot)]);
    plumesColdSlice   	= plumesCold(:,:,[UMdepthIdx(DepthIdxCold),MMdepthIdx(DepthIdxCold),LMdepthIdx(DepthIdxCold)]);
    
    %define point sources for plumes
    numPlumesHotSlice  = zeros(1,size(depthLevelSlicesPL,2));
    numPlumesColdSlice = zeros(1,size(depthLevelSlicesPL,2));
    for id=1:length(numPlumesHotSlice)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        plumesHotSlice(:,:,id)      = bwmorph(plumesHotSlice(:,:,id),'shrink',Inf);
        plumesColdSlice(:,:,id)     = bwmorph(plumesColdSlice(:,:,id),'shrink',Inf);
        try
            numPlumesHotSlice(id)    	= sum(plumesHotSlice(:,:,id),'all');
            numPlumesColdSlice(id)    	= sum(plumesColdSlice(:,:,id),'all');
        catch
            numPlumesHotSlice(id)    	= sum(sum(plumesHotSlice(:,:,id)));
            numPlumesColdSlice(id)    	= sum(sum(plumesColdSlice(:,:,id)));
        end
    end
    
    %create hot plume data arrays
    maxNumPlumesHotPerSlice         = max(numPlumesHotSlice(:));
    dummy               = zeros(size(plumesHotSlice,1),size(plumesHotSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhPlumeHotSlices    = zeros(size(depthLevelSlicesPL,2),maxNumPlumesHotPerSlice)*NaN; %[slices,plumes]
    xPlumeHotSlices     = vhPlumeHotSlices;
    yPlumeHotSlices     = vhPlumeHotSlices;
    zPlumeHotSlices     = vhPlumeHotSlices;
    for idD=1:size(depthLevelSlicesPL,2)
        %extract horizontal position
        if max(plumesHotSlice(:,:,idD))>0
            dummyNanoms = size(dummy(plumesHotSlice(:,:,idD)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(DepthIdxHot,idD));
            xPlumeHotSlices(idD,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,idD)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(DepthIdxHot,idD));
            yPlumeHotSlices(idD,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,idD)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(DepthIdxHot,idD));
            zPlumeHotSlices(idD,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,idD)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(DepthIdxHot,idD));
            vhPlumeHotSlices(idD,1:dummyNanoms)	= vhFullSlice(plumesHotSlice(:,:,idD)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    %create cold plume data arrays
    maxNumPlumesColdPerSlice         = max(numPlumesColdSlice(:));
    dummy               = zeros(size(plumesColdSlice,1),size(plumesColdSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhPlumeColdSlices   = zeros(size(depthLevelSlicesPL,2),maxNumPlumesColdPerSlice)*NaN; %[slices,plumes]
    xPlumeColdSlices    = vhPlumeColdSlices;
    yPlumeColdSlices    = vhPlumeColdSlices;
    zPlumeColdSlices    = vhPlumeColdSlices;
    for idD=1:size(depthLevelSlicesPL,2)
        %extract horizontal position
        if max(plumesColdSlice(:,:,idD))>0
            dummyNanoms = size(dummy(plumesColdSlice(:,:,idD)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(DepthIdxCold,idD));
            xPlumeColdSlices(idD,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,idD)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(DepthIdxCold,idD));
            yPlumeColdSlices(idD,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,idD)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(DepthIdxCold,idD));
            zPlumeColdSlices(idD,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,idD)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(DepthIdxCold,idD));
            vhPlumeColdSlices(idD,1:dummyNanoms)= vhFullSlice(plumesColdSlice(:,:,idD)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    %CONTINENTs
    %derive horizontal continent slices
    depthLevelSlicesCONT 	= DomainDepthMax - depthLevelCONTslices.*DomainDepth;  %plotting values
    
    depthLevelSlicesIdx = zeros(size(depthLevelSlicesCONT)); %array indices
    for id=1:length(depthLevelSlicesCONT)
        [~,depthLevelSlicesIdx(1,id)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelSlicesCONT(id)));  %index of closest value to mean LAB depth further up in the plate
    end
    
    contSlice             = CONT3D(:,:,depthLevelSlicesIdx);
    
    %define point sources for continents
    numCONTSlice = zeros(1,length(depthLevelSlicesCONT));
    for id=1:length(depthLevelSlicesCONT)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        contSlice(:,:,id)	= bwmorph(contSlice(:,:,id),'shrink',Inf);
        try
            numCONTSlice(id)  	= sum(contSlice(:,:,id),'all');
        catch
            numCONTSlice(id)  	= sum(sum(contSlice(:,:,id)));
        end
    end
    %create continent data arrays
    maxNumCONTsPerSlice= max(numCONTSlice(:));
    dummy           	= zeros(size(contSlice,1),size(contSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhCONTsSlices       = zeros(length(depthLevelSlicesCONT),maxNumCONTsPerSlice)*NaN; %[slices,continents]
    xCONTsSlices        = vhCONTsSlices;
    yCONTsSlices        = vhCONTsSlices;
    zCONTsSlices        = vhCONTsSlices;
    for id=1:length(depthLevelSlicesCONT)
        %extract horizontal position
        if max(contSlice(:,:,id))>0
            dummyNanoms = size(dummy(contSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xCONTsSlices(id,1:dummyNanoms)  	= dummy(contSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yCONTsSlices(id,1:dummyNanoms)      = dummy(contSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zCONTsSlices(id,1:dummyNanoms)      = dummy(contSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhCONTsSlices(id,1:dummyNanoms)     = vhFullSlice(contSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms

    %LLSVPs
    %derive horizontal llsvp slices
    depthLevelSlicesLLSVP 	= DomainDepthMax - depthLevelLLSVPslices.*DomainDepth;  %plotting values
    
    depthLevelSlicesIdx = zeros(size(depthLevelSlicesLLSVP)); %array indices
    for id=1:length(depthLevelSlicesLLSVP)
        [~,depthLevelSlicesIdx(1,id)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelSlicesLLSVP(id)));  %index of closest value to mean LAB depth further up in the plate
    end
    
    llsvpsSlice             = LLSVP3D(:,:,depthLevelSlicesIdx);
    
    %define point sources for llsvps
    numLLSVPSlice = zeros(1,length(depthLevelSlicesLLSVP));
    for id=1:length(depthLevelSlicesLLSVP)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        llsvpsSlice(:,:,id) 	= bwmorph(llsvpsSlice(:,:,id),'shrink',Inf);
        try
            numLLSVPSlice(id)    	= sum(llsvpsSlice(:,:,id),'all');
        catch
            numLLSVPSlice(id)    	= sum(sum(llsvpsSlice(:,:,id)));
        end
    end
    
    %create llsvp data arrays
    maxNumLLSVPsPerSlice= max(numLLSVPSlice(:));
    dummy           	= zeros(size(llsvpsSlice,1),size(llsvpsSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhLLSVPsSlices      = zeros(length(depthLevelSlicesLLSVP),maxNumLLSVPsPerSlice)*NaN; %[slices,llsvps]
    xLLSVPsSlices       = vhLLSVPsSlices;
    yLLSVPsSlices       = vhLLSVPsSlices;
    zLLSVPsSlices       = vhLLSVPsSlices;
    for id=1:length(depthLevelSlicesLLSVP)
        %extract horizontal position
        if max(llsvpsSlice(:,:,id))>0
            dummyNanoms = size(dummy(llsvpsSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xLLSVPsSlices(id,1:dummyNanoms)  	= dummy(llsvpsSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yLLSVPsSlices(id,1:dummyNanoms) 	= dummy(llsvpsSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zLLSVPsSlices(id,1:dummyNanoms) 	= dummy(llsvpsSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhLLSVPsSlices(id,1:dummyNanoms)	= vhFullSlice(llsvpsSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    %!AG! here output diagnostics can be added 
    %output diagnostics
    %Upper Mantle
    MANTLE.plumeHotNumberUM     = numPlumesHotSlice(1);
    MANTLE.plumeHotVHmaxUM      = max(vhPlumeHotSlices(1,:));
    MANTLE.plumeHotVHminUM    	= min(vhPlumeHotSlices(1,:));
    MANTLE.plumeHotVHmeanUM   	= nanmean(vhPlumeHotSlices(1,:));
    MANTLE.plumeColdNumberUM  	= numPlumesColdSlice(1);
    MANTLE.plumeColdVHmaxUM   	= max(vhPlumeColdSlices(1,:));
    MANTLE.plumeColdVHminUM    	= min(vhPlumeColdSlices(1,:));
    MANTLE.plumeColdVHmeanUM   	= nanmean(vhPlumeColdSlices(1,:));
    %Mid Mantle
    MANTLE.plumeHotNumberMM     = numPlumesHotSlice(2);
    MANTLE.plumeHotVHmaxMM      = max(vhPlumeHotSlices(2,:));
    MANTLE.plumeHotVHminMM    	= min(vhPlumeHotSlices(2,:));
    MANTLE.plumeHotVHmeanMM   	= nanmean(vhPlumeHotSlices(2,:));
    MANTLE.plumeColdNumberMM  	= numPlumesColdSlice(2);
    MANTLE.plumeColdVHmaxMM    	= max(vhPlumeColdSlices(2,:));
    MANTLE.plumeColdVHminMM    	= min(vhPlumeColdSlices(2,:));
    MANTLE.plumeColdVHmeanMM   	= nanmean(vhPlumeColdSlices(2,:));
    %Lower Mantle
    MANTLE.plumeHotNumberLM     = numPlumesHotSlice(3);
    MANTLE.plumeHotVHmaxLM      = max(vhPlumeHotSlices(3,:));
    MANTLE.plumeHotVHminLM    	= min(vhPlumeHotSlices(3,:));
    MANTLE.plumeHotVHmeanLM   	= nanmean(vhPlumeHotSlices(3,:));
    MANTLE.plumeColdNumberLM  	= numPlumesColdSlice(3);
    MANTLE.plumeColdVHmaxLM    	= max(vhPlumeColdSlices(3,:));
    MANTLE.plumeColdVHminLM    	= min(vhPlumeColdSlices(3,:));
    MANTLE.plumeColdVHmeanLM   	= nanmean(vhPlumeColdSlices(3,:));
    
    MANTLE.contVHmaxUM          = max(vhCONTsSlices(1,:));
    MANTLE.contVHminUM          = min(vhCONTsSlices(1,:));
    MANTLE.contVHmeanUM         = nanmean(vhCONTsSlices(1,:));
    
    MANTLE.llsvpVHmaxLM         = max(vhLLSVPsSlices(1,:));
    MANTLE.llsvpVHminLM         = min(vhLLSVPsSlices(1,:));
    MANTLE.llsvpVHmeanLM        = nanmean(vhLLSVPsSlices(1,:));
    
else
    warning('This has not been implemented yet for yinyang!')
end

%% DETAILS OF FOUND ANOMALIES II
% number of hot and cold active plumes
CChot2      	= bwconncomp(plumesHot);
CCcold2         = bwconncomp(plumesCold);
if strcmp(GRID.Type,'yinyang')
    CChot2_yang         	= bwconncomp(plumesHot_yang);
    CCcold2_yang         	= bwconncomp(plumesCold_yang);
    MANTLE.numHotPlumes     = CChot2.NumObjects+CChot2_yang.NumObjects;
    MANTLE.numColdPlumes    = CCcold2.NumObjects+CCcold2_yang.NumObjects;
else
    MANTLE.numHotPlumes     = CChot2.NumObjects;
    MANTLE.numColdPlumes    = CCcold2.NumObjects;
end

%check here also for wrap-around boundaries.....................


% total volume of up- and downwelling
if strcmp(GRID.Type,'yinyang')
    MANTLE.UpwellingVolume      = sum(sum(GRID.cellVolume(upWelling(:,:,2:end)==1)))+...
        sum(sum(GRID.cellVolume(upWelling_yang(:,:,2:end)==1)));     %[nd] or [m^2]or[m^3]
    MANTLE.DownwellingVolume 	= sum(sum(GRID.cellVolume(downWelling(:,:,1:end-1)==1)))+...
        sum(sum(GRID.cellVolume(downWelling_yang(:,:,1:end-1)==1))); %[nd] or [m^2]or[m^3]
else
    %!AG! in 2D spherical annulus geometry: [m^3] because of 3D scaling!
    MANTLE.UpwellingVolume      = sum(sum(GRID.cellVolume(upWelling(:,:,2:end)==1)));     %[nd] or [m^2]or[m^3]
    MANTLE.DownwellingVolume 	= sum(sum(GRID.cellVolume(downWelling(:,:,1:end-1)==1))); %[nd] or [m^2]or[m^3]
end
% percentage of total volume of up- and downwelling versus total volume model domain
MANTLE.UpwellingVolPerc         = 100 *MANTLE.UpwellingVolume /(sum(GRID.cellVolume(:))*nb);   %in [% of total volume]
MANTLE.DownwellingVolPerc  	= 100 *MANTLE.DownwellingVolume /(sum(GRID.cellVolume(:))*nb); %in [% of total volume]
% area of upwellings at certain depth..................

% active up- and downwelling
MANTLE.activeUpwelling                  = activeUpwelling;
MANTLE.activeDownwelling                = activeDownwelling;
MANTLE.activeUpwelling(~plumesHot)      = 0;
MANTLE.activeDownwelling(~plumesCold)   = 0;
CCactU                                  = bwconncomp(MANTLE.activeUpwelling);
CCactD                                  = bwconncomp(MANTLE.activeDownwelling);
MANTLE.numActUpwelling                  = CCactU.NumObjects;
MANTLE.numActDownwelling                = CCactD.NumObjects;
MANTLE.ActUpwellingVolume             	= sum(sum(GRID.cellVolume(MANTLE.activeUpwelling(:,:,2:end)==1)));      %[nd] or [m^2]or[m^3]
MANTLE.ActDownwellingVolume             = sum(sum(GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end-1)==1)));  %[nd] or [m^2]or[m^3]
MANTLE.ActUpwellingVolPerc              = 100 *MANTLE.ActUpwellingVolume /MANTLE.UpwellingVolume;               %in [% of total upwelling]
MANTLE.ActDownwellingVolPerc            = 100 *MANTLE.ActDownwellingVolume /MANTLE.DownwellingVolume;           %in [% of total downwelling]

% passive up- and downwelling
MANTLE.passiveUpwelling              	= passiveUpwelling;
MANTLE.passiveDownwelling             	= passiveDownwelling;
MANTLE.passiveUpwelling(plumesHot==1)   = 0;
MANTLE.passiveDownwelling(plumesCold==1)= 0;
CCpassU                                 = bwconncomp(MANTLE.passiveUpwelling);
CCpassD                                 = bwconncomp(MANTLE.passiveDownwelling);
MANTLE.numPassUpwelling               	= CCpassU.NumObjects;
MANTLE.numPassDownwelling               = CCpassD.NumObjects;
MANTLE.PassUpwellingVolume             	= sum(sum(GRID.cellVolume(MANTLE.passiveUpwelling(:,:,2:end)==1)));      %[nd] or [m^2]or[m^3]
MANTLE.PassDownwellingVolume            = sum(sum(GRID.cellVolume(MANTLE.passiveDownwelling(:,:,1:end-1)==1)));  %[nd] or [m^2]or[m^3]
MANTLE.PassUpwellingVolPerc             = 100 *MANTLE.PassUpwellingVolume /MANTLE.UpwellingVolume;               %in [% of total upwelling]
MANTLE.PassDownwellingVolPerc           = 100 *MANTLE.PassDownwellingVolume /MANTLE.DownwellingVolume;           %in [% of total downwelling]

% continents
if exist('CentroidCONT','var')
    MANTLE.contNumberUM                    = size(CentroidCONT,1);
    MANTLE.contLocationX                   = CentroidCONT(:,1);
    MANTLE.contLocationZ                   = CentroidCONT(:,2);
end

% llsvps
if exist('CentroidLLSVP','var')
    MANTLE.llsvpNumberLM                    = size(CentroidLLSVP,1);
    MANTLE.llsvpLocationX                   = CentroidLLSVP(:,1);
    MANTLE.llsvpLocationZ                   = CentroidLLSVP(:,2);
    %!AG!
    if strcmp(GRID.Type,'yinyang')
        MANTLE.llsvpVolume      = sum(sum(GRID.cellVolume(LLSVP3D(:,:,2:end)==1)))+...
                                  sum(sum(GRID.cellVolume(LLSVP3D_yang(:,:,2:end)==1)));     %[nd] or [m^2]or[m^3]
        MANTLE.llsvpCMBPerc         = 100 * (sum(LLSVP3D(:,1,1)==1)/(size(LLSVP3D,1)*size(LLSVP3D,2))+...
                                      (sum(LLSVP3D_yang(:,1,1)==1)/(size(LLSVP3D_yang,1)*size(LLSVP3D_yang,2))))
    else
        %!AG! in 2D spherical annulus geometry: [m^3] because of 3D scaling!
        MANTLE.llsvpVolume      = sum(sum(GRID.cellVolume(LLSVP3D(:,:,2:end)==1)));     %[nd] or [m^2]or[m^3]
        MANTLE.llsvpCMBPerc     = 100 * sum(LLSVP3D(:,1,1)==1)/(size(LLSVP3D,1)*size(LLSVP3D,2))
    end
    % percentage of total volume of up- and downwelling versus total volume model domain
    MANTLE.llsvpVolPerc         = 100 *MANTLE.llsvpVolume /(sum(GRID.cellVolume(:))*nb);   %in [% of total volume]
end


%% DETAILS OF FOUND ANOMALIES III (!AG!)
if MANTLE.numActUpwelling>0 %plume(s) present
    %Temperature
    dummyTUp                      = T_3D(MANTLE.activeUpwelling(:,:,1:end)==1);
    MANTLE.ActUpwellingT(1)       = min(dummyTUp);      % min
    MANTLE.ActUpwellingT(2)       = max(dummyTUp);      % max
    MANTLE.ActUpwellingT(3)       = sum(dummyTUp.*GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1))/MANTLE.ActUpwellingVolume;     % mean
    %Residual temperature
    dummyTresUp                   = Tres3D(MANTLE.activeUpwelling(:,:,1:end)==1);
    MANTLE.ActUpwellingTres(1)    = min(dummyTresUp);   % min
    MANTLE.ActUpwellingTres(2)    = max(dummyTresUp);   % max
    MANTLE.ActUpwellingTres(3)    = sum(dummyTresUp.*GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1))/MANTLE.ActUpwellingVolume;     % mean
    %Radial velocity
    dummyVzUp                     = VZ_3D(MANTLE.activeUpwelling(:,:,1:end)==1)*SETUP.Vscale;
    MANTLE.ActUpwellingVz(1)      = min(dummyVzUp);     % min
    MANTLE.ActUpwellingVz(2)      = max(dummyVzUp);     % max
    MANTLE.ActUpwellingVz(3)      = sum(dummyVzUp.*GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1))/MANTLE.ActUpwellingVolume;     % mean
    %Viscosity:
    dummyEtaUp                    = ETA_3D(MANTLE.activeUpwelling(:,:,1:end)==1);
    MANTLE.ActUpwellingEta(1)     = min(dummyEtaUp);    % min
    MANTLE.ActUpwellingEta(2)     = max(dummyEtaUp);    % max
    MANTLE.ActUpwellingEta(3)     = sum(dummyEtaUp.*GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1))/MANTLE.ActUpwellingVolume;     % mean
    %Viscous strain
    if exist('VSTR_3D','var')
        dummyVStrUp                   = VSTR_3D(MANTLE.activeUpwelling(:,:,1:end)==1);
        MANTLE.ActUpwellingVstr(1)    = min(dummyVStrUp);   % min
        MANTLE.ActUpwellingVstr(2)    = max(dummyVStrUp);   % max
        MANTLE.ActUpwellingVstr(3)    = sum(dummyVStrUp.*GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1))/MANTLE.ActUpwellingVolume;     % mean
        clearvars dummyVStrUp
    end
    clearvars dummyTUp dummyTresUp dummyVzUp dummyEtaUp 
end

if MANTLE.numActDownwelling>0 %slab(s) present
    %Temperature
    dummyTDown                    = T_3D(MANTLE.activeDownwelling(:,:,1:end)==1);
    MANTLE.ActDownwellingT(1)     = min(dummyTDown);    % min
    MANTLE.ActDownwellingT(2)     = max(dummyTDown);    % max
    MANTLE.ActDownwellingT(3)     = sum(dummyTDown.*GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1))/MANTLE.ActDownwellingVolume;     % mean
    %Residual temperature
    dummyTresDown                 = Tres3D(MANTLE.activeDownwelling(:,:,1:end)==1);
    MANTLE.ActDownwellingTres(1)  = min(dummyTresDown); % min
    MANTLE.ActDownwellingTres(2)  = max(dummyTresDown); % max
    MANTLE.ActDownwellingTres(3)  = sum(dummyTresDown.*GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1))/MANTLE.ActDownwellingVolume;   % mean
    %Radial velocity:
    dummyVzDown                   = VZ_3D(MANTLE.activeDownwelling(:,:,1:end)==1)*SETUP.Vscale;
    MANTLE.ActDownwellingVz(1)    = min(dummyVzDown);   % min
    MANTLE.ActDownwellingVz(2)    = max(dummyVzDown);   % max
    MANTLE.ActDownwellingVz(3)    = sum(dummyVzDown.*GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1))/MANTLE.ActDownwellingVolume;     % mean
    %Viscosity:
    dummyEtaDown                  = ETA_3D(MANTLE.activeDownwelling(:,:,1:end)==1);
    MANTLE.ActDownwellingEta(1)   = min(dummyEtaDown);  % min
    MANTLE.ActDownwellingEta(2)   = max(dummyEtaDown);  % max
    MANTLE.ActDownwellingEta(3)   = sum(dummyEtaDown.*GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1))/MANTLE.ActDownwellingVolume;     % mean
    if exist('VSTR_3D','var') %Viscous strain
        dummyVStrDown                 = VSTR_3D(MANTLE.activeDownwelling(:,:,1:end)==1);
        MANTLE.ActDownwellingVstr(1)  = min(dummyVStrDown); % min
        MANTLE.ActDownwellingVstr(2)  = max(dummyVStrDown); % max
        MANTLE.ActDownwellingVstr(3)  = sum(dummyVStrDown.*GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1))/MANTLE.ActDownwellingVolume;     % mean
        clearvars dummyVStrDown
    end
    clearvars dummyTDown dummyTresDown dummyVzDown dummyEtaDown 
end

%% DISPLAY INFORMATION
disp('   Horizontal Mantle Flow')
disp('     Velocity')
disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.VHmeanUM(DepthIdxHot),3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxUM(DepthIdxHot),3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(DepthIdxHot,1),4),' ',GRID.Zdim,' depth'])
disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.VHmeanMM(DepthIdxHot),3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxMM(DepthIdxHot),3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(DepthIdxHot,2),4),' ',GRID.Zdim,' depth'])
disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.VHmeanLM(DepthIdxHot),3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxLM(DepthIdxHot),3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(DepthIdxHot,3),4),' ',GRID.Zdim,' depth'])

disp('   Mantle Upwelling')
disp('     Volume')
%!AG! in spherical2D: volume is still in m^3 
disp(['     ',STYLE.SCHAR.smallBullet,' total                 = ',num2str(MANTLE.UpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.UpwellingVolPerc,2),' %vol)'])
disp(['     ',STYLE.SCHAR.smallBullet,' active                = ',num2str(MANTLE.ActUpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.ActUpwellingVolPerc,2),' %)'])
disp(['     ',STYLE.SCHAR.smallBullet,' passive               = ',num2str(MANTLE.PassUpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.PassUpwellingVolPerc,2),' %)'])
%     if MANTLE.numHotPlumes==0
%         disp('   No Active Hot Plume');
%     elseif MANTLE.numHotPlumes==1
%         disp('   1 Active Hot Plume:');
%     else
%         disp(['   ',num2str(MANTLE.numHotPlumes),' Active Hot Plumes:']);
%     end
if MANTLE.plumeHotNumberLM>0
    disp('     Hot-plume mobility')
    if MANTLE.plumeHotNumberUM>0
        if MANTLE.plumeHotNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.plumeHotVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberUM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxHot,1),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeHotNumberMM>0
        if MANTLE.plumeHotNumberMM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.plumeHotVHminMM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanMM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxMM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberMM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxHot,2),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeHotNumberLM>0
        if MANTLE.plumeHotNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.plumeHotVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberLM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxHot,3),4),' ',GRID.Zdim,' depth)'])
    end
end
if MANTLE.numActUpwelling>0 %!AG!
    disp('     Hot-plume diagnostics')
    disp(['     ',STYLE.SCHAR.smallBullet, ' Temp.                 = ', num2str(MANTLE.ActUpwellingT(1),4),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActUpwellingT(3),4),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActUpwellingT(2),4),' ',SETUP.TDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Residual Temp.        = ', num2str(MANTLE.ActUpwellingTres(1),4),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActUpwellingTres(3),4),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActUpwellingTres(2),4),' ',SETUP.TDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Radial velocity       = ', num2str(MANTLE.ActUpwellingVz(1),2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActUpwellingVz(3),2),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActUpwellingVz(2),2),' ',SETUP.vDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Viscosity             = ', num2str(MANTLE.ActUpwellingEta(1),2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActUpwellingEta(3),2),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActUpwellingEta(2),2),' ',SETUP.etaDim]) 
    if exist('VSTR_3D','var')  
        disp(['     ',STYLE.SCHAR.smallBullet, ' Viscous strain        = ', num2str(MANTLE.ActUpwellingVstr(1),3),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActUpwellingVstr(3),3),...
              ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActUpwellingVstr(2),3)]) 
    end
end

disp('   Mantle Downwelling')
disp('     Volume')
disp(['     ',STYLE.SCHAR.smallBullet,' total                 = ',num2str(MANTLE.DownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.DownwellingVolPerc,2),' %vol)'])
disp(['     ',STYLE.SCHAR.smallBullet,' active                = ',num2str(MANTLE.ActDownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.ActDownwellingVolPerc,2),' %)'])
%     if MANTLE.numColdPlumes==0
%         disp('   No Active Cold Plume');
%     elseif MANTLE.numColdPlumes==1
%         disp('   1 Active Cold Plume:');
%     else
%         disp(['   ',num2str(MANTLE.numColdPlumes),' Active Cold Plumes:']);
%     end
disp(['     ',STYLE.SCHAR.smallBullet,' passive               = ',num2str(MANTLE.PassDownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.PassDownwellingVolPerc,2),' %)'])

if MANTLE.plumeColdNumberLM>0
    disp('     Cold-plume mobility')
    if MANTLE.plumeColdNumberUM>0
        if MANTLE.plumeColdNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.plumeColdVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberUM,1),' slab',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxCold,1),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeColdNumberMM>0
        if MANTLE.plumeColdNumberMM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.plumeColdVHminMM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanMM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxMM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberMM,1),' slab',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxCold,2),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeColdNumberLM>0
        if MANTLE.plumeColdNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.plumeColdVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberLM,1),' slab',pluralS,' at ',num2str(depthLevelSlicesPL(DepthIdxCold,3),4),' ',GRID.Zdim,' depth)'])
    end
end
if MANTLE.numActDownwelling>0 %!AG!
    disp('     Cold-plume diagnostics')
    disp(['     ',STYLE.SCHAR.smallBullet, ' Temp.                 = ', num2str(MANTLE.ActDownwellingT(1),4),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActDownwellingT(3),4),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActDownwellingT(2),4),' ',SETUP.TDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Residual Temp.        = ', num2str(MANTLE.ActDownwellingTres(1),4),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActDownwellingTres(3),4),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActDownwellingTres(2),4),' ',SETUP.TDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Radial velocity       = ', num2str(MANTLE.ActDownwellingVz(1),2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActDownwellingVz(3),2),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActDownwellingVz(2),2),' ',SETUP.vDim])
    disp(['     ',STYLE.SCHAR.smallBullet, ' Viscosity             = ', num2str(MANTLE.ActDownwellingEta(1),2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActDownwellingEta(3),2),...
          ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActDownwellingEta(2),2),' ',SETUP.etaDim]) 
    if exist('VSTR_3D','var')  
        disp(['     ',STYLE.SCHAR.smallBullet, ' Viscous strain        = ', num2str(MANTLE.ActDownwellingVstr(1),3),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.ActDownwellingVstr(3),3),...
              ' ', STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.ActDownwellingVstr(2),3)]) 
    end
end

if MANTLE.contNumberUM>0
    disp('   Continents')
    disp('     Continent mobility')
    if MANTLE.contNumberUM>0
        if MANTLE.contNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.contVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.contVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.contVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.contNumberUM,1),' Continent',pluralS,' at ',num2str(depthLevelSlicesCONT(1),4),' ',GRID.Zdim,' depth)'])
    end
end

if MANTLE.llsvpNumberLM>0
    disp('   LLSVPs')
    if MANTLE.llsvpNumberLM>0
        if MANTLE.llsvpNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.llsvpVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.llsvpVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.llsvpVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.llsvpNumberLM,1),' LLSVP',pluralS,' at ',num2str(depthLevelSlicesLLSVP(1),4),' ',GRID.Zdim,' depth)'])
        disp(['     ',STYLE.SCHAR.smallBullet,' total                 = ',num2str(MANTLE.llsvpVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.llsvpVolPerc,2),' %vol)'])
        disp(['     ',STYLE.SCHAR.smallBullet,' CMB coverage          = ',num2str(MANTLE.llsvpCMBPerc,3),' %'])
    end
end



%% PLOTTING PLUMES SEPARATELY
if MD.plotPLUMES && strcmp(GRID.Dim,'3-D')
    figure(22),clf
    
    x = 1:size(Plumes3D,1);
    y = 1:size(Plumes3D,2);
    z = 1:size(Plumes3D,3);
    [x3d,y3d,z3d] = meshgrid(y,x,z);
    
    phot = patch(isosurface(x3d,y3d,z3d,plumesHot,+0.95)); %HOT
    isonormals(x3d,y3d,z3d,plumesHot,phot)
    set(phot,'FaceColor','red','EdgeColor','none');
    view(3);
    camlight
    lighting phong  %gouraud
    
    hold on
    pcold = patch(isosurface(x3d,y3d,z3d,plumesCold,+0.95)); %COLD
    isonormals(x3d,y3d,z3d,plumesCold,pcold)
    set(pcold,'FaceColor','blue','EdgeColor','none');
    view(3);
    axis([0 max(max(max(x3d))) 0 max(max(max(y3d))) 0 max(max(max(z3d)))])
    camlight
    lighting phong  %gouraud
    
    title('after applying criterion')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    box on
    grid on
    
    figure(23),clf
    isosurface(Plumes3D,+0.95,'noshare') %HOT
    hold on
    isosurface(Plumes3D,-0.95,'noshare') %COLD
    
    title('original')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    axis([0 max(max(max(x3d))) 0 max(max(max(y3d))) 0 max(max(max(z3d)))])
    box on
    grid on
    
    figure(1)  %go back to figure(1)
end


%% FUNCTION OUTPUT
MANTLE.upWelling            = upWelling;    %1 for upwelling, 0 else  (depending on threshold)
MANTLE.downWelling          = downWelling;  %1 for downwelling, 0 else  (depending on threshold)
MANTLE.upWellingAbsolute    = upWellingAbsolute;    %1 for upwelling, 0 else
MANTLE.downWellingAbsolute 	= downWellingAbsolute;  %1 for downwelling, 0 else
MANTLE.plumesHot            = plumesHot;    %1 for hot plumes, 0 else
MANTLE.plumesCold           = plumesCold;   %1 for cold plumes, 0 else
% MANTLE.numHotPlumes               %num
% MANTLE.numColdPlumes              %num
% MANTLE.UpwellingVolume            %[nd] or [m^2]or[m^3]
% MANTLE.DownwellingVolume          %[nd] or [m^2]or[m^3]
% MANTLE.UpwellingVolPerc           %in [% of model domain]
% MANTLE.DownwellingVolPerc         %in [% of model domain]

% MANTLE.activeUpwelling            %1 for active upwelling
% MANTLE.activeDownwelling        	%1 for active downwelling
% MANTLE.numActUpwelling           	%num
% MANTLE.numActDownwelling          %num
% MANTLE.ActUpwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.ActDownwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.ActUpwellingVolPerc    	%in [% of total upwelling]
% MANTLE.ActDownwellingVolPerc     	%in [% of total downwelling]

% MANTLE.passiveUpwelling           %1 for passive upwelling
% MANTLE.passiveDownwelling        	%1 for passive downwelling
% MANTLE.numPassUpwelling           %num
% MANTLE.numPassDownwelling         %num
% MANTLE.PassUpwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.PassDownwellingVolume      %[nd] or [m^2]or[m^3]
% MANTLE.PassUpwellingVolPerc    	%in [% of total upwelling]
% MANTLE.PassDownwellingVolPerc     %in [% of total downwelling]

if strcmp(GRID.Type,'yinyang')
    %add here......
end






















