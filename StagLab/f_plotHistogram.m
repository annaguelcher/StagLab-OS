
%%                                PLOT HISTOGRAM OF FIELD 1.0
%
% 
%                                 Anna Gülcher, 09.02.2021
%
%% NOTES
% currently just a histogram at a given timestep of a given field 
% can be mantle domain - dependent! (i.e., plumes, slabs, ..)
% multiple histograms in one plot possible (for e.g. different mantle domains)
                    
function [PLOT,SAVE] = f_plotHistogram(FIELD_M,FIELD,FILE,GRID,SWITCH,PLOT,STYLE,SETUP,SAVE,MANTLE)
%% INPUT

numberBins           = PLOT.histNumberBins;
%depthlevel          = PLOT.histDepthLevel;       % set depth level to diagnose

%% DEFAULTS
if (PLOT.histFieldLog) %!AG! or the log(XX) thing?
    normalisation    = 'pdf'; %'pdf'; %but it doesn't add up to one
else
    normalisation    = 'probability';             %'count', 'probability', 'pdf', ..
end
%'probability': relative probability, sum of bars is less or equal to 1
%               v{i}=c{i}/N
%'pdf': probability density function estimate
%               v{i}=c{i}/(N*w{i})
%               w{i} is widt of bin

%% READ INPUT DATA 
if ~strcmp(PLOT.histField,'Velocity')   %Scalar Field Data
    disp('scalar field data')    
    if strcmp(PLOT.histField,'Residual temperature')
        VAR_3D                  = zeros(size(GRID.Z_3D));
        SWITCH.customFieldTask  = 'Horizontal residual';
        [VAR]                   = f_makeField(1,FILE,GRID,SETUP,SWITCH,PLOT,1,1);
        VAR_3D(:,1,:)           = VAR.var2d;
        clearvars VAR
    elseif strcmp(PLOT.histField,'Radial velocity') %!AG!i
        %velocity
        VAR_3D                      = zeros(size(GRID.Z_3D)); %!AG! maybe not needed?
        DATA.Task                   = 'ImportFieldData';
        DATA.Field2Import           = 'Velocity';
        DATA.FieldAbbreviation      = 'V';
        [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
        if strcmp(GRID.Type,'yinyang')
            VAR_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
        else %all other grid types
            VAR_3D = PLOT.VZ_3D;
        end
        if DATA.NotFound  %!AG! this DATA.NotFound doesn't work!!
            PLOT.SFVunsuccessful = true; warning(VX_3D);
            return
        end

    else
        DATA.Task                   = 'ImportFieldData';
        DATA.Field2Import           = PLOT.histField; %see for field options in SL_FieldPlot
        DATA.FieldAbbreviation      = 'VAR';
        [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
        if strcmp(GRID.Type,'yinyang')
            VAR_3D = PLOT.VAR_3Dyin; VAR_3Dyang = PLOT.VAR_3Dyang;
        else %all other grid types
            VAR_3D = PLOT.VAR_3D;
        end
        if DATA.NotFound %!AG! this DATA.NotFound doesn't work!!
            PLOT.SFVunsuccessful = true; warning(VAR_3D);
            return
        end
    end

else %Vector Field Data: only implemented for velocity
    disp('vector field data')
    if ~strcmp(PLOT.histField,'Velocity')   %!AG! make sfvfield histField???
        PLOT.SFVunsuccessful    = true;
        warning('not implemented yet: Check here!')
        return
    end
    %velocity
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Velocity';
    DATA.FieldAbbreviation      = 'V';
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        VX_3D = PLOT.VX_3Dyin; VX_3Dyang = PLOT.VX_3Dyang;
        VY_3D = PLOT.VY_3Dyin; VY_3Dyang = PLOT.VY_3Dyang;
        VZ_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
        V_3D = PLOT.V_3Dyin; V_3Dyang = PLOT.V_3Dyang;
    else %all other grid types
        VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D; V_3D = PLOT.V_3D;
    end
    if DATA.NotFound  %!AG! this DATA.NotFound doesn't work!!
        PLOT.SFVunsuccessful = true; warning(VX_3D);
        return
    end
    if strcmp(PLOT.histField,'Velocity')
        VAR_3D      = V_3D;     %absolute velocity
        disp('VAR3D = v_3D')
        if strcmp(GRID.Type,'yinyang'); VAR_3Dyang = V_3Dyang; end
    end
end

%% LOOPING OVER MANTLE DOMAINS (histType):
for iHistLocal=1:size(PLOT.histType,2) % allow for multiple histograms for different domains

    if ~strcmp(PLOT.histType(iHistLocal), 'rprof') 
        %% EXTRACT NECCESARY DATA 
        % & COMPUTE VOLUMES; WEIGHTS 
        if strcmp(PLOT.histType(iHistLocal), 'Field')
            Data(:,:,:)         = VAR_3D;
            totalVolume       	= sum(GRID.cellVolume(:));
            w                  	= GRID.cellVolume/totalVolume;

        elseif strcmp(PLOT.histType(iHistLocal), 'FieldLM') 
            depthLevel            	= 660000; %*GRID.nd2dim
            dummy                 	= abs((GRID.Z_3D(1,1,:))-depthLevel);
            [~,depthLevelIdx]       = min(dummy); %index of closest value
            depthLevel              = GRID.Z_3Dp(1,1,depthLevelIdx); %for check 
            clearvars dummy

            Data(:,:,:)         = VAR_3D(:,:,1:depthLevelIdx);
            totalVolume         = sum(sum(GRID.cellVolume(:,:,1:depthLevelIdx)));
            w                  	= GRID.cellVolume(:,:,1:depthLevelIdx)/totalVolume;

        elseif strcmp(PLOT.histType(iHistLocal), 'FieldUM') 
            depthLevel            	= 660000; %*GRID.nd2dim
            dummy                 	= abs((GRID.Z_3D(1,1,:))-depthLevel);
            [~,depthLevelIdx]       = min(dummy); %index of closest value
            depthLevel              = GRID.Z_3Dp(1,1,depthLevelIdx); %for check 
            clearvars dummy

            Data(:,:,:)         = VAR_3D(:,:,depthLevelIdx:end);
            totalVolume         = sum(sum(GRID.cellVolume(:,:,depthLevelIdx:end)));
            w                  	= GRID.cellVolume(:,:,depthLevelIdx:end)/totalVolume;    

        elseif strcmp(PLOT.histType(iHistLocal), 'FieldPlumes')
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            Data(:,:,:)         = VAR_3D(MANTLE.activeUpwelling(:,:,1:end)==1);
            totalVolume         = sum(sum(GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1)));
            w                   = GRID.cellVolume(MANTLE.activeUpwelling(:,:,1:end)==1)/totalVolume; 

        elseif strcmp(PLOT.histType(iHistLocal), 'FieldSlabs')
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            Data(:,:,:)         = VAR_3D(MANTLE.activeDownwelling(:,:,1:end)==1);
            totalVolume         = sum(sum(GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1)));
            w                   = GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end)==1)/totalVolume;

        elseif strcmp(PLOT.histType(iHistLocal), 'FieldLLSVPs')
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            if ~(MANTLE.llsvpNumberLM>0); error('no LLSVPs detected - cannot make histogram of LLSVP field'); end
            Data(:,:,:)         = VAR_3D(MANTLE.llsvp(:,:,1:end)==1);
            totalVolume         = sum(sum(GRID.cellVolume(MANTLE.llsvp(:,:,1:end)==1)));
            w                   = GRID.cellVolume(MANTLE.llsvp(:,:,1:end)==1)/totalVolume;
  
        elseif strcmp(PLOT.histType(iHistLocal), 'FieldLLSVPs_ScaledMantle') %!AG! this doesn't work yet as weighting factor w is not taken into account into the histogram making... 
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            if ~(MANTLE.llsvpNumberLM>0); error('no LLSVPs detected - cannot make histogram of LLSVP field'); end
            Data(:,:,:)         = VAR_3D(MANTLE.llsvp(:,:,1:end)==1);
            totalVolume         = sum(sum(GRID.cellVolume(:)))
            w                   = GRID.cellVolume(MANTLE.llsvp(:,:,1:end)==1)/totalVolume;
            % !AG! but w is for now only used in the mean
            % TO DO: scale according to volume 
        elseif strcmp(PLOT.histType(iHistLocal), 'nonActive')
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            nonActiveMantle     = ones(size(VAR_3D))-MANTLE.activeDownwelling-MANTLE.activeUpwelling;
            Data(:,:,:)         = VAR_3D(nonActiveMantle==1);
            totalVolume         = sum(sum(GRID.cellVolume(nonActiveMantle==1)));
            w                   = GRID.cellVolume(nonActiveMantle==1)/totalVolume;

            clearvars nonActiveMantle

        elseif strcmp(PLOT.histType(iHistLocal), 'nonActiveLM')  
            if ~exist('MANTLE','var'); error('Mantle Diagnostics was not successful!'); end
            depthLevel            	= 660000; %*GRID.nd2dim
            dummy                 	= abs((GRID.Z_3D(1,1,:))-depthLevel);
            [~,depthLevelIdx]       = min(dummy); %index of closest value
            depthLevel              = GRID.Z_3Dp(1,1,depthLevelIdx); %for check 

            nonActiveMantle     = ones(size(VAR_3D))-MANTLE.activeDownwelling-MANTLE.activeUpwelling;
            Data(:,:,:)         = VAR_3D(nonActiveMantle(:,:,1:depthLevelIdx)==1);
            totalVolume         = sum(sum(GRID.cellVolume(nonActiveMantle(:,:,1:depthLevelIdx)==1)));
            w                   = GRID.cellVolume(nonActiveMantle(:,:,1:depthLevelIdx)==1)/totalVolume;

            clearvars dummy nonActiveMantle

        end

        %% LOG field check
        if PLOT.histFieldLog    
            Data(:,:,:) = log10(Data);
        end
        % if velocity -> extract RMS? 

        %% DIMENSINALISATION
        %!AG! is this needed?
        if ~strcmp(PLOT.histField,'Residual temperature') 
            fieldDim            = FIELD_M{strcmp(FIELD_M(:,2),PLOT.histField),4};
            varScale            = FIELD_M{strcmp(FIELD_M(:,2),PLOT.histField),5};
            Data                = Data .*varScale;
        end

        %% SQUEEZE DATA 
        %just works for 2D data 
        Data = squeeze(Data);
        w    = squeeze(w);

        %% DIAGNOSE DATA 
        DataMinValue               = min(Data(:));
        DataMaxValue               = max(Data(:));
        %weighting of cells in claculating mean, etc. 
        DataMeanValue              = sum(sum(Data.*w)/sum(w)); 
        DataMedianValue            = f_weightedMedian(Data,w);
        DataStandardDeviation      = std(Data(:),w(:));

        %WEIGHTING of histogram --> still to do! 

        %% PLOTTING
        plotAdditionLegend      = logical(1);
        if strcmp(PLOT.histType(iHistLocal), 'FieldPlumes')
            faceColor               = [0.87 0.39 0.24];
            faceAlpha               = 0.5;
            edgeColor               = [1 1 1];
            if size(PLOT.histType,2)==1 
                LineMeanColor       = [1 1 1];
            else
                LineMeanColor       = [0.94 0.82 0.78];
            end
            AreaColor               = [0.94 0.82 0.78];
            AreaAlpha               = 0.5;      
        elseif strcmp(PLOT.histType(iHistLocal), 'FieldSlabs')
            faceColor               = [0.15 0.33 0.66];
            faceAlpha               = 0.5;
            edgeColor               = [1 1 1];
            if size(PLOT.histType,2)==1 
                LineMeanColor       = [1 1 1];
            else
                LineMeanColor       = [0.76 0.84 0.96];
            end
            AreaColor               = [0.76 0.84 0.96];
            AreaAlpha               = 0.5;    
        elseif strcmp(PLOT.histType(iHistLocal), 'FieldLLSVPs') || strcmp(PLOT.histType(iHistLocal), 'FieldLLSVPs_ScaledMantle') 
            faceColor               = [0.5 0.0 0.5];
            faceAlpha               = 0.5;
            edgeColor               = [1 1 1];
            if size(PLOT.histType,2)==1
                LineMeanColor       = [1 1 1];
            else
                LineMeanColor       = [0.76 0.84 0.96];
            end
            AreaColor               = [0.76 0.84 0.96];
            AreaAlpha               = 0.5;  
        else
            faceColor               = [0.2 0.2 0.2];
            faceAlpha               = 0.5;
            edgeColor               = [1 1 1];
            if size(PLOT.histType,2)==1 
                LineMeanColor       = [1 1 1];
            else
                LineMeanColor       = [0.8 0.8 0.8];
            end
            AreaColor               = [0.9 0.9 0.9];
            AreaAlpha               = 0.6;
        end

        if strcmp(STYLE.ColorMode,'dark')
            faceColor           = 1-faceColor;
            edgeColor           = 1-edgeColor;
            LineMeanColor       = 1-LineMeanColor;
            AreaColor           = 1-AreaColor;
        end


        if isnumeric(PLOT.histXmin) && isnumeric(PLOT.histXmax)
            xMin    = PLOT.histXmin;
            xMax    = PLOT.histXmax;
        elseif isnumeric(PLOT.histXmin)
            xMin    = PLOT.histXmin;
            xMax    = DataMaxValue;
        elseif isnumeric(PLOT.histXmax)
            xMin    = DataMinValue;
            xMax    = PLOT.histXmax;
        else
            xMin    = DataMinValue;
            xMax    = DataMaxValue;
        end
        xSpan       = xMax-xMin;
        dBin        = xSpan/(numberBins);

        if PLOT.histShowOutliers
            upperLim = DataMaxValue;
        else
            upperLim = xMax;
        end

        hHist = histogram(Data,[xMin:dBin:xMax-dBin,upperLim],...
            'FaceColor',faceColor,'FaceAlpha',faceAlpha,'EdgeColor',edgeColor,'EdgeAlpha',0.8,'Normalization',normalisation);
        hold on
        %constant axis
        if SWITCH.Colorbar(1,2)
            if isnumeric(PLOT.histYmax)
                ylim(gca,[PLOT.histYmin,PLOT.histYmax])
            else
                ylim(gca,[PLOT.histYmin,FIELD.cColorbarMAX])
            end
        end
        if isnumeric(PLOT.histXmin) || isnumeric(PLOT.histXmax)
            xlim(gca,[xMin,xMax])
        end
        if PLOT.histLogY %!AG!
            set(gca, 'YScale', 'log')
        end

        %plot additions
        if PLOT.histAdditions
            hold on
            yAxisLim    = ylim;
            if size(PLOT.histType,2)==1  %only plott stdv when 1 histogram is shown
                %standard deviation
                hStdDev = fill([DataMeanValue-DataStandardDeviation,DataMeanValue-DataStandardDeviation,...
                    DataMeanValue+DataStandardDeviation,DataMeanValue+DataStandardDeviation],...
                    [yAxisLim,fliplr(yAxisLim)],AreaColor,'EdgeColor','none','FaceAlpha',AreaAlpha);
                hold on
            end
            %mean value
            try
                LineMeanColor = [LineMeanColor,0.8]; %transparency
                hMean = plot([DataMeanValue,DataMeanValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7);
            catch
                LineMeanColor(1,4) = [];
                hMean = plot([DataMeanValue,DataMeanValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7);
            end
            %median value
            hMedian = plot([DataMedianValue,DataMedianValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7,'LineStyle',':');

            %put histogram to the top %!AG! somehow only works without X axis inlogscale
            uistack(hHist,'top') 
            if plotAdditionLegend
                %legend
                if size(PLOT.histType,2)==1
                    legend([hStdDev,hMean,hMedian],{'Std. Deviation','Mean','Median'});
                else
                    legend([hMean,hMedian],{'Mean','Median'});
                end
            end
        end

        %% ANNOTIATIONS
        if PLOT.histFieldLog
            DESVARIA.xlabelName         = ['log(', PLOT.histField,')'];
        else
            DESVARIA.xlabelName         = PLOT.histField;
        end
        if strcmp(PLOT.histField,'Residual temperature')
            DESVARIA.xlabelDim              = 'K';
        else
            DESVARIA.xlabelDim              = fieldDim;
        end
        if strcmp(normalisation,'probability') || strcmp(normalisation,'pdf')
            DESVARIA.zlabelName         = 'Probability';
            DESVARIA.zlabelDim          = '';
        elseif strcmp(normalisation,'counts')
            DESVARIA.zlabelName         = 'Frequency';
            DESVARIA.zlabelDim          = '';
        else
            DESVARIA.zlabelName         = '';
            DESVARIA.zlabelDim          = '';
            warning('Normalisation not recognised: Check here!')
        end

        DESVARIA.Task   = 'create annotation strings';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

        title(PLOT.titleStringCurrent)
        xlabel(PLOT.xlabel)
        ylabel(PLOT.ylabel)
        set(gca,'YGrid','on')
        set(gca,'TickDir',SWITCH.TickDirection)

        if (iHistLocal==1)
            DESVARIA.Task   = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        end

        %% FUNCTION OUTPUT
        PLOT.histTotalVolume             = [PLOT.histTotalVolume,    totalVolume    ];
        PLOT.histDataMinValue            = [PLOT.histDataMinValue,   DataMinValue   ];
        PLOT.histDataMaxValue            = [PLOT.histDataMaxValue,   DataMaxValue   ];
        PLOT.histDataMeanValue           = [PLOT.histDataMeanValue,  DataMeanValue  ];
        PLOT.histDataMedianValue         = [PLOT.histDataMedianValue,DataMedianValue];
        PLOT.histDataStandardDeviation   = [PLOT.histDataStandardDeviation,DataStandardDeviation];
        %PLOT.histDepthLevel              = depthLevel;       %update depth level in plotting dimension

    else 
        %% RPROF histogram
        %!AG!: TO DO

    end
    
    clearvars Data totalVolume DataMinValue DataMaxValue DataMeanValue DataMedianValue DataStandardDeviation
end %end loop over histogram fields 
%end

end


%% INTERIOR FUNCTIONS
%
function [wMed] = f_weightedMedian(D,W)
% ----------------------------------------------------------------------
% Function for calculating the weighted median 
% Sven Haase, altered by AG 2021
%
% For n numbers x_1,...,x_n with positive weights w_1,...,w_n, 
% (sum of all weights equal to one) the weighted median is defined as
% the element x_k, such that:
%           --                        --
%           )   w_i  <= 1/2   and     )   w_i <= 1/2
%           --                        --
%        x_i < x_k                 x_i > x_k
%
%
% Input:    D ... matrix of observed values
%           W ... matrix of weights, W = ( w_ij )
% Output:   wMed ... weighted median                   
% ----------------------------------------------------------------------
if nargin ~= 2
    error('weightedMedian:wrongNumberOfArguments', ...
      'Wrong number of arguments.');
end
if size(D) ~= size(W)
    error('weightedMedian:wrongMatrixDimension', ...
      'The dimensions of the input-matrices must match.');
end
%normalize the weights, such that: sum ( w_ij ) = 1
% (sum of all weights equal to one)
WSum = sum(W(:));
W = W / WSum;
%check dimensions of data and weight arrays
NdimD = ndims(D);
NdimW = ndims(W);
if ~(NdimD==NdimW) %!AG!
    error('dimensions of data and weight arrays not equal')
end
if NdimD>2 %!AG!
    % 3D arrays: should be 2D
    % remove dimension of size 1
    D = squeeze(D);
    W = squeeze(W);
    disp([ 'size of data: '])
    size(D)
    disp([ 'size of w   : '])
    size(W)
end

% (line by line) transformation of the input-matrices to line-vectors
d = reshape(D',1,[]);   
w = reshape(W',1,[]);  
% sort the vectors
A = [d' w'];
ASort = sortrows(A,1);
dSort = ASort(:,1)';
wSort = ASort(:,2)';
sumVec = [];    % vector for cumulative sums of the weights
for i = 1:length(wSort)
    sumVec(i) = sum(wSort(1:i));
end
wMed = [];      
j = 0;         
while isempty(wMed)
    j = j + 1;
    if sumVec(j) >= 0.5
        wMed = dSort(j);    % value of the weighted median
    end
end
% final test to exclude errors in calculation
if ( sum(wSort(1:j-1)) > 0.5 ) & ( sum(wSort(j+1:length(wSort))) > 0.5 )
     error('weightedMedian:unknownError', ...
      'The weighted median could not be calculated.');
end

end


