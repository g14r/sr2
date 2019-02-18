function sr2_training_makeTargetFiles(varargin)
%%       sr2_training_makeTargetFiles(varargin)
%
% rt: reaction time task
% st: sequence training task
% rs: random sequence task
%
% inputs: subject number vector, and session/day number vector
% outputs: target files for each subject: rt task (2 blocks), st task (10 blocks), rs task (2 blocks)
%
% example call:
%               sr2_training_makeTargetFiles([1:2],[1:4]);
%
% --
% gariani@uwo.ca - 2017.10.12

%% make target files for all subjects (varargin{1}) and sessions/days (varargin{2})
for s=varargin{1}
    for d=varargin{2}
        fprintf(1,'\nsubj: %d   session/day: %d\n',s,d);
        sr2_rt_target(s,d);
        sr2_strs_target(s,d);
    end
end
end

function[varargout]=sr2_rt_target(s,d)
%%      [varargout]=sr2_rt_target(s,d)
%
% rt: reaction time task
%
% inputs: subject number s, and session/day number d
% outputs: target files for rt task (2 blocks)
%
% standalone usage:
%        [varargout]=sr2_rt_target(1,1);
%
% --
% gariani@uwo.ca - 2017.10.12

%% define target folder
targetFolder='../../../../robotcode/projects/SequenceRepetition/sr2/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist

%% experimental details per session
seq={
    [1 8 7 0 9]; %1
    [1 9 0 8 7]; %2
    [2 9 8 6 0]; %3
    [2 0 6 8 9]; %4
    [3 6 7 9 0]; %5
    [3 0 9 6 7]; %6
    [4 6 0 7 8]; %7
    [4 7 8 0 6]; %8
    [5 7 6 9 8]; %9
    [5 8 9 7 6] %10
    };

nSeq=size(seq,1); %number of sequences
prepTime=[
    200; 250; 300; 350; 400;
    450; 500; 550; 600; 650
    ];
nPrepTime=numel(prepTime);
n_rt_blocks=2;
respWindow=100;
nTrialsPerBlock=nSeq*5;
nTrials_rt=nTrialsPerBlock*n_rt_blocks;
fixedTrialDuration=0; %1|0 use fixed trial duration (fMRI study), or make it variable depending on subjects' movement time (behavioral study)

%% define dataframe structures for rt task
Drt.repType=ones(nTrials_rt,1);
Drt.seqNum=floor((((1:nTrials_rt)-1)/(nTrials_rt/nPrepTime)+1)');
Drt.repNum=ones(nTrials_rt,1);
Drt.prepTime=repmat(prepTime,nTrials_rt/nPrepTime,1);
Drt.respWindow=repmat(respWindow,nTrials_rt,1);
Drt.dummy=zeros(nTrials_rt,1);
Drt.train=zeros(nTrials_rt,1);
Drt.rtt=ones(nTrials_rt,1);
Drt.day=ones(nTrials_rt,1)*d;

%% perform randomization
randIdx=randperm(nTrials_rt);%randomize indices for seqNum order
Drt=getrow(Drt,randIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
Drt.BN=floor((((1:nTrials_rt)-1)/(nTrialsPerBlock)+1)'+n_rt_blocks*(d-1)); %add BN info

%% split Drt structure into different block structures
tStart=0; %start time for the first trial
trialDur=6000; %approximate duration of each trial (in ms)
itiDur=500; %ITI duration (in ms)
filename=struct(); %preallocate structure (one per block)

Brt=cell(n_rt_blocks,1); %preallocate B structure (one per block)
for b=1:n_rt_blocks
    Brt{b}=getrow(Drt,Drt.BN==b+n_rt_blocks*(d-1));
    
    %add variable number of warmup (dummy) trials (repType=1, prepTime=3) at the beginning of each block
    nWarmupTrials=5;
    if nWarmupTrials>0
        %add warmup trials at the beginning of each block
        warmup.repType=ones(nWarmupTrials,1);
        warmup.seqNum=randperm(nSeq,nWarmupTrials)';
        warmup.repNum=ones(nWarmupTrials,1);
        warmup.BN=ones(nWarmupTrials,1)*b+n_rt_blocks*(d-1);
        warmup.prepTime=ones(nWarmupTrials,1)*prepTime(end);
        warmup.respWindow=ones(nWarmupTrials,1)*respWindow;
        warmup.dummy=ones(nWarmupTrials,1);
        warmup.train=zeros(nWarmupTrials,1);
        warmup.rtt=ones(nWarmupTrials,1);
        warmup.day=ones(nWarmupTrials,1)*d;
        %merge structs
        Brt{b}=insertrow(Brt{b},1,warmup); %beginning of block
    end
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Brt{b}.seqNum)*(trialDur+itiDur);
        Brt{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %rtart time for each trial
    else
        Brt{b}.tStart=ones(numel(Brt{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Brt{b}.iti=ones(numel(Brt{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Brt{b}.feedback=ones(numel(Brt{b}.seqNum),1); %give feedback (yes/no)
    Brt{b}.hand=ones(numel(Brt{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Brt{b}.seqNum)
        Brt{b}.press1(t,1)=seq{Brt{b}.seqNum(t)}(1);
        Brt{b}.press2(t,1)=seq{Brt{b}.seqNum(t)}(2);
        Brt{b}.press3(t,1)=seq{Brt{b}.seqNum(t)}(3);
        Brt{b}.press4(t,1)=seq{Brt{b}.seqNum(t)}(4);
        Brt{b}.press5(t,1)=seq{Brt{b}.seqNum(t)}(5);
        Brt{b}.cueP{t,1}=num2str(seq{Brt{b}.seqNum(t)},'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Brt{b}=rmfield(Brt{b},extraFields);
    
    %% save structure Brt{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_rt_s%02d_b%02d.tgt',s,b+n_rt_blocks*(d-1)));
    dsave(filename.fn{b,1},Brt{b})
end

%% return output data structure D
varargout{1}=filename;
varargout{2}=Brt;
end


function[varargout]=sr2_strs_target(s,d,varargin)
%%      [varargout]=sr2_strs_target(s,d,varargin)
%
% st: sequence training task
% rs: random sequence task
%
% inputs: subject number s, and session/day number d
% outputs: target files for both st and rs tasks (10 and 2 blocks, respectively)
%
% standalone usage:
%        [varargout]=sr2_strs_target(1,1);
%
%
% --
% gariani@uwo.ca - 2017.10.12

%% define target folder
targetFolder='../../../../robotcode/projects/SequenceRepetition/sr2/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist

%% experimental details per session
% find all possible sequences with each finger only once
allseq=perms([1,2,3,4,5]);

% remove "runs" (e.g. 1,2,3, ... 5,4,3, ...)
allseq_noRun=zeros(1,size(allseq,2));
rowCount=1;
for r=1:size(allseq,1)
    isRun=zeros(1,size(allseq,2)-2);
    for c=1:size(allseq,2)-2
        if allseq(r,c)==(allseq(r,c+1)-1) && allseq(r,c)==(allseq(r,c+2)-2)
            isRun(1,c)=1;
        elseif allseq(r,c)==(allseq(r,c+1)+1) && allseq(r,c)==(allseq(r,c+2)+2)
            isRun(1,c)=1;
        else
            isRun(1,c)=0;
        end
    end
    if all(isRun==0)
        allseq_noRun(rowCount,:)=allseq(r,:);
        rowCount=rowCount+1;
    else
        continue
    end
end

% pick train set
train_set={
    [1 3 2 4 5]; %1 A1 11
    [1 4 5 2 3]; %2 A2 12
    [2 4 3 5 1]; %3 B1 21
    [2 5 1 4 3]; %4 B2 22
    [3 1 2 5 4]; %5 C1 31
    [3 5 4 2 1]; %6 C2 32
    [4 1 5 3 2]; %7 D1 41
    [4 2 3 1 5]; %8 D2 42
    [5 2 1 3 4]; %9 E1 51
    [5 3 4 1 2] %10 E2 52
    };
% train_set={
%     [1 3 2 4 5]; %1 
%     [2 5 3 4 1]; %2 
%     [3 5 2 1 4]; %3 
%     [5 4 2 3 1]; %4 
%     };

% check first-order transitions in train set
trans_mat_train=zeros(5,5);
for r=1:size(train_set,1)
    for c=1:size(train_set{1},2)-1
        trans_mat_train(train_set{r}(c),train_set{r}(c+1))=trans_mat_train(train_set{r}(c),train_set{r}(c+1))+1;
    end
end
%trans_mat_train_pct=(trans_mat_train/sum(sum(trans_mat_train)))*100;

% find left-out test set
train_idx=false(size(allseq_noRun,1),1);
for i=1:size(allseq_noRun,1)
    for j=1:size(train_set,1)
        if isequal(allseq_noRun(i,:),train_set{j})
            train_idx(i,1)=true;
        end
    end
end
test_idx=~train_idx;
test_set=allseq_noRun(test_idx,:);

% randomize and remove 2 sequences
test_set=test_set(randperm(size(test_set,1)),:);
test_set=test_set(1:end-2,:);
%test_set=test_set(1:end,:);

% check first-order transitions in test set
trans_mat_test=zeros(5,5);
for r=1:size(test_set,1)
    for c=1:size(test_set,2)-1
        trans_mat_test(test_set(r,c),test_set(r,c+1))=trans_mat_test(test_set(r,c),test_set(r,c+1))+1;
    end
end
%trans_mat_test_pct=(trans_mat_test/sum(sum(trans_mat_test)))*100;

nSeq_train=size(train_set,1); %number of sequences
nSeq_test=size(test_set,1); %number of sequences
prepTime=[400; 800; 1600; 2400];
nPrepTime=numel(prepTime);
n_st_blocks=10;
n_rs_blocks=2;
respWindow=100;
nTrialsPerBlock=nSeq_train*nPrepTime;
nTrials_st=nTrialsPerBlock*n_st_blocks;
nTrials_rs=nSeq_test;
fixedTrialDuration=0; %1|0 use fixed trial duration (fMRI study), or make it variable depending on subjects' movement time (behavioral study)

%% define dataframe structures for st and rs tasks
Dst.repType=ones(nTrials_st,1);
Dst.seqNum=floor((((1:nTrials_st)-1)/(nTrialsPerBlock)+1)');
Dst.repNum=ones(nTrials_st,1);
Dst.prepTime=repmat(prepTime,nSeq_train*n_st_blocks,1);
Dst.respWindow=repmat(respWindow,nTrials_st,1);
Dst.dummy=zeros(nTrials_st,1);
Dst.train=ones(nTrials_st,1);
Dst.rtt=zeros(nTrials_st,1);
Dst.day=ones(nTrials_st,1)*d;

Drs.repType=ones(nTrials_rs,1);
Drs.seqNum=(1:nTrials_rs)';
Drs.repNum=ones(nTrials_rs,1);
Drs.prepTime=repmat(prepTime,nTrials_rs/nPrepTime,1);
Drs.respWindow=repmat(respWindow,nTrials_rs,1);
Drs.dummy=zeros(nTrials_rs,1);
Drs.train=zeros(nTrials_rs,1);
Drs.rtt=zeros(nTrials_rs,1);
Drs.day=ones(nTrials_rs,1)*d;

%% perform randomization
randIdx=randperm(nTrials_st);%randomize indices for seqNum order
Dst=getrow(Dst,randIdx); %implement randomization
Dst.BN=floor((((1:nTrials_st)-1)/(nTrialsPerBlock)+1)')+n_st_blocks*(d-1); %add BN info

randIdx=randperm(nTrials_rs);%randomize indices for seqNum order
Drs=getrow(Drs,randIdx); %implement randomization
Drs.BN=floor((((1:nTrials_rs)-1)/(nTrialsPerBlock)+1)')+n_rs_blocks*(d-1); %add BN info

%% split Dst and Drs structures into different block structures Bst and Brs
tStart=0; %start time for the first trial
trialDur=6000; %approximate duration of each trial (in ms)
itiDur=500; %ITI duration (in ms)
filename=struct(); %preallocate structure (one per block)

Bst=cell(n_st_blocks,1); %preallocate B structure (one per block)
for b=1:n_st_blocks
    Bst{b}=getrow(Dst,Dst.BN==b+n_st_blocks*(d-1));
    
    %add variable number of warmup (dummy) trials (repType=1, prepTime=3) at the beginning of each block
    nWarmupTrials=5;
    if nWarmupTrials>0
        %add warmup trials at the beginning of each block
        warmup.repType=ones(nWarmupTrials,1);
        warmup.seqNum=randperm(nSeq_train,nWarmupTrials)';
        warmup.repNum=ones(nWarmupTrials,1);
        warmup.BN=ones(nWarmupTrials,1)*b+n_st_blocks*(d-1);
        warmup.prepTime=ones(nWarmupTrials,1)*prepTime(end);
        warmup.respWindow=ones(nWarmupTrials,1)*respWindow;
        warmup.dummy=ones(nWarmupTrials,1);
        warmup.train=zeros(nWarmupTrials,1);
        warmup.rtt=zeros(nWarmupTrials,1);
        warmup.day=ones(nWarmupTrials,1)*d;
        %merge structs
        Bst{b}=insertrow(Bst{b},1,warmup); %beginning of block
    end
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Bst{b}.seqNum)*(trialDur+itiDur);
        Bst{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Bst{b}.tStart=ones(numel(Bst{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Bst{b}.iti=ones(numel(Bst{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Bst{b}.feedback=ones(numel(Bst{b}.seqNum),1); %give feedback (yes/no)
    Bst{b}.hand=ones(numel(Bst{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Bst{b}.seqNum)
        Bst{b}.press1(t,1)=train_set{Bst{b}.seqNum(t)}(1);
        Bst{b}.press2(t,1)=train_set{Bst{b}.seqNum(t)}(2);
        Bst{b}.press3(t,1)=train_set{Bst{b}.seqNum(t)}(3);
        Bst{b}.press4(t,1)=train_set{Bst{b}.seqNum(t)}(4);
        Bst{b}.press5(t,1)=train_set{Bst{b}.seqNum(t)}(5);
        Bst{b}.cueP{t,1}=num2str(train_set{Bst{b}.seqNum(t)},'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Bst{b}=rmfield(Bst{b},extraFields);
    
    %% save structure Bst{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_st_s%02d_b%02d.tgt',s,b+n_st_blocks*(d-1)));
    dsave(filename.fn{b,1},Bst{b})
end

Brs=cell(n_rs_blocks,1); %preallocate B structure (one per block)
for b=1:n_rs_blocks
    Brs{b}=getrow(Drs,Drs.BN==b+n_rs_blocks*(d-1));
    
    %add variable number of warmup (dummy) trials (repType=1, prepTime=3) at the beginning of each block
    nWarmupTrials=5;
    if nWarmupTrials>0
        %add warmup trials at the beginning of each block
        warmup.repType=ones(nWarmupTrials,1);
        warmup.seqNum=randperm(nSeq_test,nWarmupTrials)';
        warmup.repNum=ones(nWarmupTrials,1);
        warmup.BN=ones(nWarmupTrials,1)*b+n_rs_blocks*(d-1);
        warmup.prepTime=ones(nWarmupTrials,1)*prepTime(end);
        warmup.respWindow=ones(nWarmupTrials,1)*respWindow;
        warmup.dummy=ones(nWarmupTrials,1);
        warmup.train=zeros(nWarmupTrials,1);
        warmup.rtt=zeros(nWarmupTrials,1);
        warmup.day=ones(nWarmupTrials,1)*d;
        %merge structs
        Brs{b}=insertrow(Brs{b},1,warmup); %beginning of block
    end
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Brs{b}.seqNum)*(trialDur+itiDur);
        Brs{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Brs{b}.tStart=ones(numel(Brs{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Brs{b}.iti=ones(numel(Brs{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Brs{b}.feedback=ones(numel(Brs{b}.seqNum),1); %give feedback (yes/no)
    Brs{b}.hand=ones(numel(Brs{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Brs{b}.seqNum)
        Brs{b}.press1(t,1)=test_set(Brs{b}.seqNum(t),1);
        Brs{b}.press2(t,1)=test_set(Brs{b}.seqNum(t),2);
        Brs{b}.press3(t,1)=test_set(Brs{b}.seqNum(t),3);
        Brs{b}.press4(t,1)=test_set(Brs{b}.seqNum(t),4);
        Brs{b}.press5(t,1)=test_set(Brs{b}.seqNum(t),5);
        Brs{b}.cueP{t,1}=num2str(test_set(Brs{b}.seqNum(t),:),'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Brs{b}=rmfield(Brs{b},extraFields);
    
    %% save structure Brs{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_rs_s%02d_b%02d.tgt',s,b+n_rs_blocks*(d-1)));
    dsave(filename.fn{b,1},Brs{b})
end

%% return output data structure D
varargout{1}=filename;
varargout{2}=Bst;
varargout{3}=Brs;
end