function sr2_makeTargetFiles(varargin)
%% sr1_makeTargetFiles(varagin)
% example call:
%               sr2_makeTargetFiles([1:2],[1]);
%               sr2_makeTargetFiles([1:2],[1],'rw');
%               sr2_makeTargetFiles([1:2],[1],'pt');
%
% --
% gariani@uwo.ca - 2017.09.05
%% make target files for all subjects (varargin{1}) and runs (varargin{2}) per subject specified by vector inputs in varargin
for s=varargin{1}
    for r=varargin{2}
        fprintf(1,'\nsubj: %d   run: %d\n',s,r);
        if numel(varargin)==2; [fn,~]=sr2_target(s,r); else; [fn,~]=sr2_target(s,r,varargin(3)); end
        while isempty(fn)
            fprintf(1,'Trying again...\n\nsubj: %d   run: %d\n',s,r)
            if numel(varargin)==2; [fn,~]=sr2_target(s,r); else; [fn,~]=sr2_target(s,r,varargin(3)); end
        end
    end
end
end

function[varargout]=sr2_target(s,r,varargin)
%% [varargout]=sr2_target(s,r)
% Creates .tgt files (one per block, four blocks) for subject s and run r
% This function generates four separate .tgt files (one per block, within
% run r) for each time that it is called
% inputs: subject number (s), run number (r), whether it is a pilot ('rw','pt') or not
% outputs: saved filename (fn), block structure (B)
%
% standalone usage:
%        [fn,B,T]=sr2_target(1,1);
%        [fn,B,T]=sr2_target(1,1,'rw');
%        [fn,B,T]=sr2_target(1,1,'pt');
%
% for multiple subjects, call the fucntion:
%        sr2_makeTargetFiles(1,1);
%
% --
% gariani@uwo.ca - 2017.09.19

%% define target folder
cd ../../../..; targetFolder='robotcode/projects/SequenceRepetition/sr2/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist
cd matlab/project/SequenceRepetition/sr2 %go back to original matlab code folder

%% experimental details per run
fixedTrialDuration=0; %1|0 use fixed trial duration (fMRI study), or make it variable depending on subjects' movement time (behavioral study)
% allSeq={...   %equals to perms([1,2,3,5])
%     [1 2 3 5];
%     [1 2 5 3];
%     [1 3 5 2];
%     [1 3 2 5];
%     [1 5 2 3];
%     [1 5 3 2];
%     [2 1 3 5];
%     [2 1 5 3];
%     [2 3 1 5];
%     [2 3 5 1];
%     [2 5 1 3];
%     [2 5 3 1];
%     [3 1 2 5];
%     [3 1 5 2];
%     [3 2 1 5];
%     [3 2 5 1];
%     [3 5 1 2];
%     [3 5 2 1];
%     [5 1 2 3];
%     [5 1 3 2];
%     [5 2 3 1];
%     [5 2 1 3];
%     [5 3 1 2];
%     [5 3 2 1]};
seq={... %chosen set of four-digit sequences (each finger twice in each position, no more than 2 neighboring fingers in a row)
    [1 3 2 5];   %1
    [1 5 2 3];   %2
    [2 3 1 5];   %3
    [2 5 1 3];   %4
    [3 1 5 2];   %5
    [3 2 5 1];   %6
    [5 1 3 2];   %7
    [5 2 3 1]};  %8
nSeq=numel(seq); %number of sequences
nMaxRep=3; %max number of same sequence repetitions
respWindowLevels=100;
nRespWindow=numel(respWindowLevels);
prepTimeLevels=[4000,4000,4000]; %[250,1000,4000]; %manipulate preparation time (3 levels, in ms)
nPrepTime=numel(prepTimeLevels);
p=0.25; q=1-p; prb=zeros(1,nMaxRep);
for i=1:nMaxRep
    prb(i)=p*q^(i-1); %probability of sequence repetitions (geometric distribution with p=0.25)
end
pct=prb/sum(prb); %proportion of trials per repetition type
nRepPerRepType=round(pct*10); %number of repetitions in a run per sequence repetition type
nTrepType=1:nMaxRep; %number of trials per repetition type
nTrialsPerRun=sum(nRepPerRepType.*nTrepType)*nSeq*nPrepTime; %number of trials in each run

%% define dataframe structures for data (D) and trial (T)
D.repType=[];
for j=1:numel(nTrepType); D.repType=[D.repType;ones(nRepPerRepType(j),1)*j]; end %compute how many trials for each repetition type, per sequence
D.repType=repmat(D.repType,nSeq*nPrepTime,1); %extend count to full run (all sequences)
D.seqNum=floor((((1:(sum(nRepPerRepType)*nSeq*nPrepTime))-1)/((sum(nRepPerRepType)*nSeq*nPrepTime)/nSeq))+1)'; %ordered sequence type vector

%% perform randomization
seqOrderIdx=randperm(numel(D.seqNum));%randomize indices for seqNum order
tic %start timer
for k=1:numel(seqOrderIdx)-1 %for each index
    while D.seqNum(seqOrderIdx(k))==D.seqNum(seqOrderIdx(k+1)) %check whether next index points to the same seqNum (seqNum repetition)
        %if it does, keep shuffling the rest of the seqOrderIdx vector until it doesn't anymore
        if toc>3 %if more than 3s have passed (i.e. infinite loop due to impossible shuffling) display warning, exit function with empty outputs, and try again
            warning('entering infinite loop due to impossible shuffling! No files saved: returning control to the invoking function.')
            varargout{1}=[]; %empty output 1
            varargout{2}=[]; %empty output 2
            return
        end
        temp=seqOrderIdx(k+1:end);
        idx=randperm(numel(temp));
        temp=temp(idx);
        seqOrderIdx(k+1:end)=temp;
    end
end
D=getrow(D,seqOrderIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D

%% fill in data structure T
T.repType=[];
T.repTypeCount=[];
T.seqNum=[];
T.repNum=[];
for l=1:numel(D.repType) %for each repType, expand fields of T
    if D.repType(l)==1
        T.repType=[T.repType;ones(1,1)*D.repType(l)]; %zeros and NaNs correspond to announce trials at the beginning of each new sequence type
        T.repTypeCount=[T.repTypeCount;ones(1,1)*l];
        T.seqNum=[T.seqNum;ones(1,1)*D.seqNum(l)];
        T.repNum=[T.repNum;1];
    elseif D.repType(l)==2
        T.repType=[T.repType;ones(2,1)*D.repType(l)];
        T.repTypeCount=[T.repTypeCount;ones(2,1)*l];
        T.seqNum=[T.seqNum;ones(2,1)*D.seqNum(l)];
        T.repNum=[T.repNum;[1;2]];
    elseif D.repType(l)==3
        T.repType=[T.repType;ones(3,1)*D.repType(l)];
        T.repTypeCount=[T.repTypeCount;ones(3,1)*l];
        T.seqNum=[T.seqNum;ones(3,1)*D.seqNum(l)];
        T.repNum=[T.repNum;[1;2;3]];
    elseif D.repType(l)==4
        T.repType=[T.repType;ones(4,1)*D.repType(l)];
        T.repTypeCount=[T.repTypeCount;ones(4,1)*l];
        T.seqNum=[T.seqNum;ones(4,1)*D.seqNum(l)];
        T.repNum=[T.repNum;[1;2;3;4]];
    elseif D.repType(l)==5
        T.repType=[T.repType;ones(5,1)*D.repType(l)];
        T.repTypeCount=[T.repTypeCount;ones(5,1)*l];
        T.seqNum=[T.seqNum;ones(5,1)*D.seqNum(l)];
        T.repNum=[T.repNum;[1;2;3;4;5]];
    else
        error('Exceeded max number of same sequence repetitions! Check nMaxRep')
    end
end

%% is it a pilot or not? What kind of pilot in case?
if ~isempty(varargin); varargin=varargin{1}; end 
if strcmp(varargin,'rw'); respWindowPilot=1; prepTimePilot=0;
    respWindowLevels=25:25:200;
    nRespWindow=numel(respWindowLevels);
    prepTimeLevels=4000;
    nPrepTime=numel(prepTimeLevels);
elseif strcmp(varargin,'pt'); respWindowPilot=0; prepTimePilot=1;
    respWindowLevels=100;
    nRespWindow=numel(respWindowLevels);
    prepTimeLevels=50:50:400;
    nPrepTime=numel(prepTimeLevels);
else
    respWindowPilot=0; prepTimePilot=0;
end
T.prepTime=(ceil((1:nTrialsPerRun)/(nTrialsPerRun/nPrepTime)))'; for i=1:nPrepTime; T.prepTime(T.prepTime==i)=prepTimeLevels(i); end %ordered prepTime vector
T.prepTime=T.prepTime(randperm(nTrialsPerRun)); %randomized prepTime vector
T.respWindow=(ceil((1:nTrialsPerRun)/(nTrialsPerRun/nRespWindow)))'; for i=1:nRespWindow; T.respWindow(T.respWindow==i)=respWindowLevels(i); end %ordered respWindow vector
T.respWindow=T.respWindow(randperm(nTrialsPerRun)); %randomized respWindow vector
T.TN=(1:nTrialsPerRun)'; %add info on how many trials per run

%% break down each structure run (R) into 8 block structures (B) of about 48 trials each without chopping off repType
%  to do so, add buffer trials at the beginning and end of each block (50 trials fixed)
%initialize/preallocate count variables
approxNtrialsPerBlock=50;
nBlocks=ceil(nTrialsPerRun/approxNtrialsPerBlock);
trialCount=1;
trialCountVec=zeros(1,nBlocks);
blockNum=1;
startBlock=1;
blockIndices=cell(1,nBlocks);
R=struct();
for trialNum=1:numel(T.TN) %for each trial, check whether it's end of block or not
    if trialCount>(approxNtrialsPerBlock-nMaxRep) && T.repNum(trialNum)==1 && blockNum<nBlocks %every about 47(min)-49(max) trials, end of block (for blocks 1-7)
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
tStart=0; %start time for the first trial
trialDur=6000; %approximate duration of each trial (in ms)
itiDur=500; %ITI duration (in ms)
B=cell(nBlocks,1); %preallocate B structure (one per block)
filename=struct(); %preallocate structure (one per block)
cd ../../../..;
for b=1:nBlocks
    B{b}=getrow(R,R.BN==b);
    
    %add variable number of warmup (dummy) trials (repType=1, prepTime=3) at the beginning of each block
    nWarmupTrials=approxNtrialsPerBlock-numel(B{b}.seqNum);
    if nWarmupTrials>0
        %add warmup trials at the beginning of each block
        warmup.repType=ones(nWarmupTrials,1);
        warmup.repTypeCount=zeros(nWarmupTrials,1);
        warmup.seqNum=randperm(nSeq,nWarmupTrials)'; while warmup.seqNum(end)==B{b}.seqNum(1); warmup.seqNum=randperm(nSeq,nWarmupTrials)'; end %make sure there is no seqNum repetition between the last warmup trial and the first actual trial
        warmup.repNum=ones(nWarmupTrials,1);
        warmup.TN=zeros(nWarmupTrials,1);
        warmup.TNblock=zeros(nWarmupTrials,1);
        warmup.BN=ones(nWarmupTrials,1)*b;
        warmup.dummy=ones(nWarmupTrials,1);
        warmup.prepTime=ones(nWarmupTrials,1)*prepTimeLevels(end);
        warmup.respWindow=ones(nWarmupTrials,1)*respWindowLevels(end);
        %merge structs
        B{b}=insertrow(B{b},1,warmup); %beginning of block
    end
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(B{b}.seqNum)*(trialDur+itiDur);
        B{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        B{b}.tStart=ones(numel(B{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
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
    extraFields={'BN','TN','TNblock','repTypeCount'};
    B{b}=rmfield(B{b},extraFields);
    
    %% save structure B{b} as a target file (.tgt)
    if prepTimePilot==1
        filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2pt_s%02d_b%02d.tgt',s,8*(r-1)+b));
    elseif respWindowPilot==1
        filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2rw_s%02d_b%02d.tgt',s,8*(r-1)+b));
    else
        filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_s%02d_b%02d.tgt',s,8*(r-1)+b));
    end
    dsave(filename.fn{b,1},B{b})
end
cd matlab/project/SequenceRepetition/sr2 %go back to original matlab code folder

%% sanity check!
% [...
%     B{1}.BN,B{1}.TNrun,B{1}.TN,B{1}.TNblock,B{1}.dummy,B{1}.tStart,B{1}.prepTime,B{1}.iti,B{1}.seqNum,B{1}.repType,B{1}.repNum,B{1}.feedback,B{1}.hand,B{1}.press1,B{1}.press2,B{1}.press3,B{1}.press4;...
%     B{2}.BN,B{2}.TNrun,B{2}.TN,B{2}.TNblock,B{2}.dummy,B{2}.tStart,B{2}.prepTime,B{2}.iti,B{2}.seqNum,B{2}.repType,B{2}.repNum,B{2}.feedback,B{2}.hand,B{2}.press1,B{2}.press2,B{2}.press3,B{2}.press4;...
%     B{3}.BN,B{3}.TNrun,B{3}.TN,B{3}.TNblock,B{3}.dummy,B{3}.tStart,B{3}.prepTime,B{3}.iti,B{3}.seqNum,B{3}.repType,B{3}.repNum,B{3}.feedback,B{3}.hand,B{3}.press1,B{3}.press2,B{3}.press3,B{3}.press4;...
%     B{4}.BN,B{4}.TNrun,B{4}.TN,B{4}.TNblock,B{4}.dummy,B{4}.tStart,B{4}.prepTime,B{4}.iti,B{4}.seqNum,B{4}.repType,B{4}.repNum,B{4}.feedback,B{4}.hand,B{4}.press1,B{4}.press2,B{4}.press3,B{4}.press4]

% [B{1}.TNrun,B{1}.TN,B{1}.TNblock,B{1}.dummy,B{1}.tStart,B{1}.prepTime,B{1}.iti,B{1}.seqNum,B{1}.repType,B{1}.repNum,B{1}.feedback,B{1}.hand,B{1}.press1,B{1}.press2,B{1}.press3,B{1}.press4]
% [B{2}.TNrun,B{2}.TN,B{2}.TNblock,B{2}.dummy,B{2}.tStart,B{2}.prepTime,B{2}.iti,B{2}.seqNum,B{2}.repType,B{2}.repNum,B{2}.feedback,B{2}.hand,B{2}.press1,B{2}.press2,B{2}.press3,B{2}.press4]
% [B{3}.TNrun,B{3}.TN,B{3}.TNblock,B{3}.dummy,B{3}.tStart,B{3}.prepTime,B{3}.iti,B{3}.seqNum,B{3}.repType,B{3}.repNum,B{3}.feedback,B{3}.hand,B{3}.press1,B{3}.press2,B{3}.press3,B{3}.press4]
% [B{4}.TNrun,B{4}.TN,B{4}.TNblock,B{4}.dummy,B{4}.tStart,B{4}.prepTime,B{4}.iti,B{4}.seqNum,B{4}.repType,B{4}.repNum,B{4}.feedback,B{4}.hand,B{4}.press1,B{4}.press2,B{4}.press3,B{4}.press4]

%% return output data structure D
varargout{1}=filename;
varargout{2}=B;
varargout{3}=T;
end