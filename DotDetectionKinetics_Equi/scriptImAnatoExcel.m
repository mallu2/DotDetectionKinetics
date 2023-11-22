analysisDir = '/crex/proj/uppstore2018129/elflab/Projects/LacIdynamicsInVivo/EXP-23-CA4324/';
addpath(fullfile(analysisDir,'2309_code'));
addpath(fullfile(analysisDir,'2309_code','ImAnalysis'));
ImAnalysis_setup(fullfile(analysisDir,'2309_code','ImAnalysis'))
% figDir = fullfile(analysisDir,'figures');
% mkdir(figDir)
%% Run
dataDir = '/domus/h1/malinlu/elflab/Projects/LacIdynamicsInVivo/EXP-23-CA4325/';

% specifiy the run
run = 'run10_2';
% the parameter file 'even' or 'odd'
side = 'odd';
parameterfile = [analysisDir 'ParameterFile_CA4325_' side '.txt'];
 
sourceDir = [dataDir '/the' run];
outputDir = fullfile(analysisDir,'output',run,side);
mkdir(outputDir)

processCells(sourceDir, outputDir, parameterfile);

% Get the number of dots per cell
expInfoObj = loadExpInfo(outputDir);
fluoChanName = expInfoObj.getChannelNames('fluo');
fluoChanName = fluoChanName{1};


% initiate array for times
PositionTime = zeros(1,length(expInfoObj.positions));

%get the position list
posList = expInfoObj.getPositionList();
dotInfo = [];
parfor p=1:numel(expInfoObj.positions)
    mCells = expInfoObj.getMCells(p);
    dotdata = load(fullfile(outputDir,posList{p},['LSQdotCoordinates_' fluoChanName '.mat']),'allMAPpar','totalDotCoord0');
    dotList = getParticleListForMalin(mCells,dotdata.allMAPpar,dotdata.totalDotCoord0);
    dotList(1).Position=repmat(p,size(dotList(1).coord,1),1);
    dotInfo = [dotInfo;  dotList];
    
end

% Export data to Excel
fields=fieldnames(dotInfo);
dl = struct();
for f=1:numel(fields)
    dl.(fields{f}) =vertcat(dotInfo.(fields{f}));
end
tbl = struct2table(dl);
writetable(tbl,fullfile(analysisDir,['dots_', run '_' side, '.xlsx']),'WriteMode','replacefile');
