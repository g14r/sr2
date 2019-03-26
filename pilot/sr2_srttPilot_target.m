function[varargout]=sr2_srttPilot_target(s,r,varargin)
%% [varargout]=sr2_srttPilot_target(s,r)
% Creates .tgt files for subject s and run r
% inputs: subject number (s), run number (r)
% outputs: saved filename (fn), block structure (B)
%
% standalone usage:
%        [fn,B,T]=sr2_srttPilot_target(1,1);
%
% --
% gariani@uwo.ca - 2017.09.21

%% define target folder
cd ../../../..; targetFolder='robotcode/projects/SequenceRepetition/sr2/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist
cd matlab/project/SequenceRepetition/sr2 %go back to original matlab code folder

%% experimental details per run
seq={...
    [1 7 6 8];   %1
    [1 8 9 7];   %2
    [2 6 9 7];   %3
    [2 9 7 8];   %4
    [3 8 6 7];   %5
    [3 7 9 8];   %6
    [5 9 8 7];   %7
    [5 6 9 8]};  %8
respWindowLevels=100;
nRespWindow=numel(respWindowLevels);
% prepTimeLevels=[100;125;150;175;
%                 200;225;250;275;
%                 300;325;350;375;
%                 400;450;500;800];
% prepTimeLevels=[
%     125;150;175;200;
%     225;250;275;300;
%     325;350;375;400;
%     500;800;1600;2400
%     ];

prepTimeLevels=[
    150;175;200;225;
    250;275;300;325;
    350;375;400;450;
    500;800;1600;2400
    ];

nTrialsPerRun=384;
one=repmat((1:2)',8,1);
two=repmat((3:4)',8,1);
three=repmat((5:6)',8,1);
five=repmat((7:8)',8,1);

%% define dataframe structures for data (D) and trial (T)
D.repType=ones(nTrialsPerRun,1);
D.repNum=ones(nTrialsPerRun,1);
D.seqNum=repmat([one;two;three;five],6,1);
D.prepTime=repmat(prepTimeLevels,24,1);
D.respWindow=(ceil((1:nTrialsPerRun)/(nTrialsPerRun/nRespWindow)))'; for i=1:nRespWindow; D.respWindow(D.respWindow==i)=respWindowLevels(i); end

%% perform randomization
seqOrderIdx=randperm(numel(D.seqNum));%randomize indices for seqNum order
T=getrow(D,seqOrderIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
T.TN=(1:nTrialsPerRun)'; %add info on how many trials per run

%% break down each structure run (R) into 8 block structures (B) of about 48 trials each without chopping off repType
approxNtrialsPerBlock=48;
nBlocks=ceil(nTrialsPerRun/approxNtrialsPerBlock);
trialCount=1;
trialCountVec=zeros(1,nBlocks);
blockNum=1;
startBlock=1;
blockIndices=cell(1,nBlocks);
R=struct();
for trialNum=1:numel(T.TN) %for each trial, check whether it's end of block or not
    if trialCount>approxNtrialsPerBlock && T.repNum(trialNum)==1 && blockNum<nBlocks %every about 47(min)-49(max) trials, end of block (for blocks 1-7)
        endBlock=trialNum-1;
        trialCountVec(blockNum)=numel(startBlock:endBlock);
        blockIndices{:,blockNum}=(startBlock:endBlock)'; %get the trial indices for beginning/end of this block
        thisBlock=getrow(T,blockIndices{:,blockNum});
        thisBlock.TNblock=(1:(numel(thisBlock.seqNum)))';
        thisBlock.BN=ones(numel(thisBlock.seqNum),1)*blockNum;
        thisBlock.dummy=zeros(numel(thisBlock.seqNum),1);
        R=addstruct(R,thisBlock,'row');
        trialCount=2; %reset trial count (considering current trialNum as 1)
        blockNum=blockNum+1; %move on to the next block
        startBlock=trialNum; %keep track of new startBlock index
    elseif trialNum==numel(T.TN) %end of run, block nBlocks
        endBlock=trialNum;
        trialCountVec(blockNum)=numel(startBlock:endBlock);
        blockIndices{:,blockNum}=(startBlock:endBlock)'; %get the trial indices for beginning/end of this block
        thisBlock=getrow(T,blockIndices{:,blockNum});
        thisBlock.TNblock=(1:(numel(thisBlock.seqNum)))';
        thisBlock.BN=ones(numel(thisBlock.seqNum),1)*blockNum;
        thisBlock.dummy=zeros(numel(thisBlock.seqNum),1);
        R=addstruct(R,thisBlock,'row');
    else %within blocks
        trialCount=trialCount+1; %move to the next trial
    end
end

%split run structure R into n different block structures B{b}
itiDur=500; %ITI duration (in ms)
B=cell(nBlocks,1); %preallocate B structure (one per block)
filename=struct(); %preallocate structure (one per block)
cd ../../../..;
for b=1:nBlocks
    B{b}=getrow(R,R.BN==b);
    
    %add exp timing info
    B{b}.tStart=ones(numel(B{b}.seqNum),1)*(-1); %use invalid cue -1 to
    B{b}.iti=ones(numel(B{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    B{b}.feedback=ones(numel(B{b}.seqNum),1); %give feedback (yes/no)
    B{b}.hand=ones(numel(B{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(B{b}.seqNum)
        B{b}.press1(t,1)=seq{B{b}.seqNum(t)}(1);
        B{b}.press2(t,1)=seq{B{b}.seqNum(t)}(2);
        B{b}.press3(t,1)=seq{B{b}.seqNum(t)}(3);
        B{b}.press4(t,1)=seq{B{b}.seqNum(t)}(4);
        B{b}.cueP{t,1}=num2str(seq{B{b}.seqNum(t)},'%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN','TN','TNblock'};
    B{b}=rmfield(B{b},extraFields);
    
    %% save structure B{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2srt_s%02d_b%02d.tgt',s,8*(r-1)+b));
    dsave(filename.fn{b,1},B{b})
end
cd matlab/project/SequenceRepetition/sr2 %go back to original matlab code folder

%% return output data structure D
varargout{1}=filename;
varargout{2}=B;
varargout{3}=T;
end