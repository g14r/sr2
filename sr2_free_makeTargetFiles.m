function sr2_free_makeTargetFiles(varargin)
%%       sr2_free_makeTargetFiles(varargin)
%
% fr: free-RT task
% fc: forced-RT task
%
% inputs: subject number vector, and session/day number vector
% outputs: target files for each subject: fr task (4 blocks), fc task (4 blocks)
%
% example call:
%               sr2_free_makeTargetFiles(1,1);
%
% --
% gariani@uwo.ca - 2018.01.09

%% make target files for all subjects (varargin{1}) and sessions/days (varargin{2})
for s=varargin{1}
    for d=varargin{2}
        fprintf(1,'\nsubj: %d   session/day: %d\n',s,d);
        sr2_fr_target(s,d);
        sr2_fc_target(s,d);
    end
end
end

function[varargout]=sr2_fr_target(s,d)
%%      [varargout]=sr2_fr_target(s,d)
%
% fr: free-RT task
%
% inputs: subject number s, and session/day number d
% outputs: target files for fr task (8 blocks)
%
% standalone usage:
%        [varargout]=sr2_fr_target(1,1);
%
% --
% gariani@uwo.ca - 2018.01.09

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
    [1 3 2 4 5]; %1 A1
    [1 4 5 2 3]; %2 A2
    [2 4 3 5 1]; %3 B1
    [2 5 1 4 3]; %4 B2
    [3 1 2 5 4]; %5 C1
    [3 5 4 2 1]; %6 C2
    [4 1 5 3 2]; %7 D1
    [4 2 3 1 5]; %8 D2
    [5 2 1 3 4]; %9 E1
    [5 3 4 1 2] %10 E2
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

% find left-out untrained set
train_idx=false(size(allseq_noRun,1),1);
for i=1:size(allseq_noRun,1)
    for j=1:size(train_set,1)
        if isequal(allseq_noRun(i,:),train_set{j})
            train_idx(i,1)=true;
        end
    end
end
untrained_idx=~train_idx;
untrained_set=allseq_noRun(untrained_idx,:);

% randomize and remove 2 sequences
untrained_set=untrained_set(randperm(size(untrained_set,1)),:);
untrained_set=untrained_set(1:end-2,:);

% check first-order transitions in untrained set
trans_mat_untrained=zeros(5,5);
for r=1:size(untrained_set,1)
    for c=1:size(untrained_set,2)-1
        trans_mat_untrained(untrained_set(r,c),untrained_set(r,c+1))=trans_mat_untrained(untrained_set(r,c),untrained_set(r,c+1))+1;
    end
end
%trans_mat_untrained_pct=(trans_mat_untrained/sum(sum(trans_mat_untrained)))*100;

%% TRAINING: BLOCKED (tr: trained vs un: untrained)
nSeq_train=size(train_set,1); %number of trained sequences
nSeq_untrained=size(untrained_set,1); %number of untrained sequences
prepTime=[0; 0; 0; 0];
nPrepTime=numel(prepTime);
n_tr_blocks=3;
n_un_blocks=1;
respWindow=100;
nTrialsPerBlock=nSeq_train*nPrepTime;
nTrials_tr=nTrialsPerBlock*n_tr_blocks;
nTrials_un=nSeq_untrained;
fixedTrialDuration=0; %1|0 use fixed trial duration (fMRI study), or make it variable depending on subjects' movement time (behavioral study)

%% define dataframe structures for tr and un blocks
Dtr.repType=ones(nTrials_tr,1);
Dtr.seqNum=floor((((1:nTrials_tr)-1)/(nTrials_tr/nSeq_train)+1)');
Dtr.repNum=ones(nTrials_tr,1);
Dtr.prepTime=repmat(prepTime,nSeq_train*n_tr_blocks,1);
Dtr.respWindow=repmat(respWindow,nTrials_tr,1);
Dtr.dummy=zeros(nTrials_tr,1);
Dtr.train=ones(nTrials_tr,1);
Dtr.rtt=zeros(nTrials_tr,1);
Dtr.free=ones(nTrials_tr,1);
Dtr.mix=zeros(nTrials_tr,1);
Dtr.day=ones(nTrials_tr,1)*d;

Dun.repType=ones(nTrials_un,1);
Dun.seqNum=(1:nTrials_un)'+nSeq_train;
Dun.repNum=ones(nTrials_un,1);
Dun.prepTime=repmat(prepTime,nTrials_un/nPrepTime,1);
Dun.respWindow=repmat(respWindow,nTrials_un,1);
Dun.dummy=zeros(nTrials_un,1);
Dun.train=zeros(nTrials_un,1);
Dun.rtt=zeros(nTrials_un,1);
Dun.free=ones(nTrials_un,1);
Dun.mix=zeros(nTrials_un,1);
Dun.day=ones(nTrials_un,1)*d;

%% perform randomization
randIdx=randperm(nTrials_tr);%randomize indices for seqNum order
Dtr=getrow(Dtr,randIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
Dtr.BN=floor((((1:nTrials_tr)-1)/(nTrialsPerBlock)+1)')+n_tr_blocks*(d-1); %add BN info

randIdx=randperm(nTrials_un);%randomize indices for seqNum order
Dun=getrow(Dun,randIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
Dun.BN=floor((((1:nTrials_un)-1)/(nTrialsPerBlock)+1)')+n_un_blocks*(d-1); %add BN info

%% split Dtr and Dun structures into different block structures Btr and Bun
tStart=0; %start time for the first trial
trialDur=6000; %approximate duration of each trial (in ms)
itiDur=500; %ITI duration (in ms)
filename=struct(); %preallocate structure (one per block)

Btr=cell(n_tr_blocks,1); %preallocate B structure (one per block)
for b=1:n_tr_blocks
    Btr{b}=getrow(Dtr,Dtr.BN==b+n_tr_blocks*(d-1));
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Btr{b}.seqNum)*(trialDur+itiDur);
        Btr{b}.ttrart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Btr{b}.tStart=ones(numel(Btr{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Btr{b}.iti=ones(numel(Btr{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Btr{b}.feedback=ones(numel(Btr{b}.seqNum),1); %give feedback (yes/no)
    Btr{b}.hand=ones(numel(Btr{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Btr{b}.seqNum)
        Btr{b}.press1(t,1)=train_set{Btr{b}.seqNum(t)}(1);
        Btr{b}.press2(t,1)=train_set{Btr{b}.seqNum(t)}(2);
        Btr{b}.press3(t,1)=train_set{Btr{b}.seqNum(t)}(3);
        Btr{b}.press4(t,1)=train_set{Btr{b}.seqNum(t)}(4);
        Btr{b}.press5(t,1)=train_set{Btr{b}.seqNum(t)}(5);
        Btr{b}.cueP{t,1}=num2str(train_set{Btr{b}.seqNum(t)},'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Btr{b}=rmfield(Btr{b},extraFields);
    
    %% save structure Btr{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fr_tr_s%02d_b%02d.tgt',s,b+n_tr_blocks*(d-1)));
    dsave(filename.fn{b,1},Btr{b})
end

Bun=cell(n_un_blocks,1); %preallocate B structure (one per block)
for b=1:n_un_blocks
    Bun{b}=getrow(Dun,Dun.BN==b+n_un_blocks*(d-1));
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Bun{b}.seqNum)*(trialDur+itiDur);
        Bun{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Bun{b}.tStart=ones(numel(Bun{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Bun{b}.iti=ones(numel(Bun{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Bun{b}.feedback=ones(numel(Bun{b}.seqNum),1); %give feedback (yes/no)
    Bun{b}.hand=ones(numel(Bun{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Bun{b}.seqNum)
        Bun{b}.press1(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,1);
        Bun{b}.press2(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,2);
        Bun{b}.press3(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,3);
        Bun{b}.press4(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,4);
        Bun{b}.press5(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,5);
        Bun{b}.cueP{t,1}=num2str(untrained_set(Bun{b}.seqNum(t)-nSeq_train,:),'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Bun{b}=rmfield(Bun{b},extraFields);
    
    %% save structure Bun{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fr_un_s%02d_b%02d.tgt',s,b+n_un_blocks*(d-1)));
    dsave(filename.fn{b,1},Bun{b})
end

%% TRAINING: MIXED (mx: trained and untrained)
n_mx_blocks=4;
n_trained_trials_per_block=30;
n_untrained_trials_per_block=10;
nTrialsPerBlock=n_trained_trials_per_block+n_untrained_trials_per_block;
n_trials_mx=nTrialsPerBlock*n_mx_blocks;

%% define dataframe structures for mx blocks
Dmx.repType=ones(n_trials_mx,1);
Dmx.seqNum=repmat((repmat((1:nSeq_train)',(nTrialsPerBlock/nSeq_train),1)),n_mx_blocks,1);
Dmx.repNum=ones(n_trials_mx,1);
Dmx.prepTime=repmat(prepTime,n_trials_mx/nPrepTime,1);
Dmx.respWindow=repmat(respWindow,n_trials_mx,1);
Dmx.dummy=zeros(n_trials_mx,1);
Dmx.train=ones(n_trials_mx,1);
Dmx.rtt=zeros(n_trials_mx,1);
Dmx.free=ones(n_trials_mx,1);
Dmx.mix=ones(n_trials_mx,1);
Dmx.day=ones(n_trials_mx,1)*d;
Dmx.BN=floor((((1:n_trials_mx)-1)/(nTrialsPerBlock)+1)')+n_mx_blocks*(d-1); %add BN info

%% split Dmx structure into different block structure Bmx
filename=struct(); %preallocate structure (one per block)
Bmx=cell(n_mx_blocks,1); %preallocate B structure (one per block)
for b=1:n_mx_blocks
    Bmx{b}=getrow(Dmx,Dmx.BN==b+n_mx_blocks*(d-1));
    
    %add untrained sequences (10 trials) and mix trained/untrained trials within block
    iStart=(nTrialsPerBlock+1+(n_untrained_trials_per_block*(b-1))); iEnd=(nTrialsPerBlock+n_untrained_trials_per_block*b);
    Bmx{b}.seqNum(1:n_untrained_trials_per_block,1)=Dun.seqNum(iStart:iEnd);
    Bmx{b}.train(1:n_untrained_trials_per_block,1)=zeros(n_untrained_trials_per_block,1);    
    randIdx=randperm(numel(Bmx{b}.seqNum));%randomize indices for seqNum order
    Bmx{b}=getrow(Bmx{b},randIdx); %implement randomization
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Bmx{b}.seqNum)*(trialDur+itiDur);
        Bmx{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Bmx{b}.tStart=ones(numel(Bmx{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Bmx{b}.iti=ones(numel(Bmx{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Bmx{b}.feedback=ones(numel(Bmx{b}.seqNum),1); %give feedback (yes/no)
    Bmx{b}.hand=ones(numel(Bmx{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Bmx{b}.seqNum)
        if Bmx{b}.seqNum(t)<=nSeq_train 
            %trained seq trials
            Bmx{b}.press1(t,1)=train_set{Bmx{b}.seqNum(t)}(1);
            Bmx{b}.press2(t,1)=train_set{Bmx{b}.seqNum(t)}(2);
            Bmx{b}.press3(t,1)=train_set{Bmx{b}.seqNum(t)}(3);
            Bmx{b}.press4(t,1)=train_set{Bmx{b}.seqNum(t)}(4);
            Bmx{b}.press5(t,1)=train_set{Bmx{b}.seqNum(t)}(5);
            Bmx{b}.cueP{t,1}=num2str(train_set{Bmx{b}.seqNum(t)},'%d%d%d%d%d');
        else
            %untrained seq trials
            Bmx{b}.press1(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,1);
            Bmx{b}.press2(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,2);
            Bmx{b}.press3(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,3);
            Bmx{b}.press4(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,4);
            Bmx{b}.press5(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,5);
            Bmx{b}.cueP{t,1}=num2str(untrained_set(Bmx{b}.seqNum(t)-nSeq_train,:),'%d%d%d%d%d');
        end
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Bmx{b}=rmfield(Bmx{b},extraFields);
    
    %% save structure Bmx{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fr_mx_s%02d_b%02d.tgt',s,b+n_mx_blocks*(d-1)));
    dsave(filename.fn{b,1},Bmx{b})
end

%% return output data structure D
varargout{1}=filename;
varargout{2}=Btr;
varargout{3}=Bun;
varargout{4}=Bmx;
end


function[varargout]=sr2_fc_target(s,d,varargin)
%%      [varargout]=sr2_fc_target(s,d,varargin)
%
% fc: forced-RT task
%
% inputs: subject number s, and session/day number d
% outputs: target files for the fc task (8 blocks)
% 4 blocks training blocked (3 trained 1 untrained, 40 trials per block)
% 4 blocks training mixed (30 trials trained 10 trials untrained per block)
%
% standalone usage:
%        [varargout]=sr2_fc_target(1,1);
%
%
% --
% gariani@uwo.ca - 2018.01.09

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
    [1 3 2 4 5]; %1 A1
    [1 4 5 2 3]; %2 A2
    [2 4 3 5 1]; %3 B1
    [2 5 1 4 3]; %4 B2
    [3 1 2 5 4]; %5 C1
    [3 5 4 2 1]; %6 C2
    [4 1 5 3 2]; %7 D1
    [4 2 3 1 5]; %8 D2
    [5 2 1 3 4]; %9 E1
    [5 3 4 1 2] %10 E2
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

% find left-out untrained set
train_idx=false(size(allseq_noRun,1),1);
for i=1:size(allseq_noRun,1)
    for j=1:size(train_set,1)
        if isequal(allseq_noRun(i,:),train_set{j})
            train_idx(i,1)=true;
        end
    end
end
untrained_idx=~train_idx;
untrained_set=allseq_noRun(untrained_idx,:);

% randomize and remove 2 sequences
untrained_set=untrained_set(randperm(size(untrained_set,1)),:);
untrained_set=untrained_set(1:end-2,:);

% check first-order transitions in untrained set
trans_mat_untrained=zeros(5,5);
for r=1:size(untrained_set,1)
    for c=1:size(untrained_set,2)-1
        trans_mat_untrained(untrained_set(r,c),untrained_set(r,c+1))=trans_mat_untrained(untrained_set(r,c),untrained_set(r,c+1))+1;
    end
end
%trans_mat_untrained_pct=(trans_mat_untrained/sum(sum(trans_mat_untrained)))*100;

%% TRAINING: BLOCKED (tr: trained vs un: untrained)
nSeq_train=size(train_set,1); %number of trained sequences
nSeq_untrained=size(untrained_set,1); %number of untrained sequences
prepTime=[400; 800; 1600; 2400];
nPrepTime=numel(prepTime);
n_tr_blocks=3;
n_un_blocks=1;
respWindow=100;
nTrialsPerBlock=nSeq_train*nPrepTime;
nTrials_tr=nTrialsPerBlock*n_tr_blocks;
nTrials_un=nSeq_untrained;
fixedTrialDuration=0; %1|0 use fixed trial duration (fMRI study), or make it variable depending on subjects' movement time (behavioral study)

%% define dataframe structures for tr and un blocks
Dtr.repType=ones(nTrials_tr,1);
Dtr.seqNum=floor((((1:nTrials_tr)-1)/(nTrials_tr/nSeq_train)+1)');
Dtr.repNum=ones(nTrials_tr,1);
Dtr.prepTime=repmat(prepTime,nSeq_train*n_tr_blocks,1);
Dtr.respWindow=repmat(respWindow,nTrials_tr,1);
Dtr.dummy=zeros(nTrials_tr,1);
Dtr.train=ones(nTrials_tr,1);
Dtr.rtt=zeros(nTrials_tr,1);
Dtr.free=zeros(nTrials_tr,1);
Dtr.mix=zeros(nTrials_tr,1);
Dtr.day=ones(nTrials_tr,1)*d;

Dun.repType=ones(nTrials_un,1);
Dun.seqNum=(1:nTrials_un)'+nSeq_train;
Dun.repNum=ones(nTrials_un,1);
Dun.prepTime=repmat(prepTime,nTrials_un/nPrepTime,1);
Dun.respWindow=repmat(respWindow,nTrials_un,1);
Dun.dummy=zeros(nTrials_un,1);
Dun.train=zeros(nTrials_un,1);
Dun.rtt=zeros(nTrials_un,1);
Dun.free=zeros(nTrials_un,1);
Dun.mix=zeros(nTrials_un,1);
Dun.day=ones(nTrials_un,1)*d;

%% perform randomization
randIdx=randperm(nTrials_tr);%randomize indices for seqNum order
Dtr=getrow(Dtr,randIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
Dtr.BN=floor((((1:nTrials_tr)-1)/(nTrialsPerBlock)+1)')+n_tr_blocks*(d-1); %add BN info

randIdx=randperm(nTrials_un);%randomize indices for seqNum order
Dun=getrow(Dun,randIdx); %using the optimized seqOrderIdx, implement randomization of both repType and seqNum vectors and store in structure D
Dun.BN=floor((((1:nTrials_un)-1)/(nTrialsPerBlock)+1)')+n_un_blocks*(d-1); %add BN info

%% split Dtr and Dun structures into different block structures Btr and Bun
tStart=0; %start time for the first trial
trialDur=6000; %approximate duration of each trial (in ms)
itiDur=500; %ITI duration (in ms)
filename=struct(); %preallocate structure (one per block)

Btr=cell(n_tr_blocks,1); %preallocate B structure (one per block)
for b=1:n_tr_blocks
    Btr{b}=getrow(Dtr,Dtr.BN==b+n_tr_blocks*(d-1));
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Btr{b}.seqNum)*(trialDur+itiDur);
        Btr{b}.ttrart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Btr{b}.tStart=ones(numel(Btr{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Btr{b}.iti=ones(numel(Btr{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Btr{b}.feedback=ones(numel(Btr{b}.seqNum),1); %give feedback (yes/no)
    Btr{b}.hand=ones(numel(Btr{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Btr{b}.seqNum)
        Btr{b}.press1(t,1)=train_set{Btr{b}.seqNum(t)}(1);
        Btr{b}.press2(t,1)=train_set{Btr{b}.seqNum(t)}(2);
        Btr{b}.press3(t,1)=train_set{Btr{b}.seqNum(t)}(3);
        Btr{b}.press4(t,1)=train_set{Btr{b}.seqNum(t)}(4);
        Btr{b}.press5(t,1)=train_set{Btr{b}.seqNum(t)}(5);
        Btr{b}.cueP{t,1}=num2str(train_set{Btr{b}.seqNum(t)},'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Btr{b}=rmfield(Btr{b},extraFields);
    
    %% save structure Btr{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fc_tr_s%02d_b%02d.tgt',s,b+n_tr_blocks*(d-1)));
    dsave(filename.fn{b,1},Btr{b})
end

Bun=cell(n_un_blocks,1); %preallocate B structure (one per block)
for b=1:n_un_blocks
    Bun{b}=getrow(Dun,Dun.BN==b+n_un_blocks*(d-1));
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Bun{b}.seqNum)*(trialDur+itiDur);
        Bun{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Bun{b}.tStart=ones(numel(Bun{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Bun{b}.iti=ones(numel(Bun{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Bun{b}.feedback=ones(numel(Bun{b}.seqNum),1); %give feedback (yes/no)
    Bun{b}.hand=ones(numel(Bun{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Bun{b}.seqNum)
        Bun{b}.press1(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,1);
        Bun{b}.press2(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,2);
        Bun{b}.press3(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,3);
        Bun{b}.press4(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,4);
        Bun{b}.press5(t,1)=untrained_set(Bun{b}.seqNum(t)-nSeq_train,5);
        Bun{b}.cueP{t,1}=num2str(untrained_set(Bun{b}.seqNum(t)-nSeq_train,:),'%d%d%d%d%d');
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Bun{b}=rmfield(Bun{b},extraFields);
    
    %% save structure Bun{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fc_un_s%02d_b%02d.tgt',s,b+n_un_blocks*(d-1)));
    dsave(filename.fn{b,1},Bun{b})
end

%% TRAINING: MIXED (mx: trained and untrained)
n_mx_blocks=4;
n_trained_trials_per_block=30;
n_untrained_trials_per_block=10;
nTrialsPerBlock=n_trained_trials_per_block+n_untrained_trials_per_block;
n_trials_mx=nTrialsPerBlock*n_mx_blocks;

%% define dataframe structures for mx blocks
Dmx.repType=ones(n_trials_mx,1);
Dmx.seqNum=repmat((repmat((1:nSeq_train)',(nTrialsPerBlock/nSeq_train),1)),n_mx_blocks,1);
Dmx.repNum=ones(n_trials_mx,1);
Dmx.prepTime=repmat(prepTime,n_trials_mx/nPrepTime,1);
Dmx.respWindow=repmat(respWindow,n_trials_mx,1);
Dmx.dummy=zeros(n_trials_mx,1);
Dmx.train=ones(n_trials_mx,1);
Dmx.rtt=zeros(n_trials_mx,1);
Dmx.free=zeros(n_trials_mx,1);
Dmx.mix=ones(n_trials_mx,1);
Dmx.day=ones(n_trials_mx,1)*d;
Dmx.BN=floor((((1:n_trials_mx)-1)/(nTrialsPerBlock)+1)')+n_mx_blocks*(d-1); %add BN info

%% split Dmx structure into block structure Bmx
filename=struct(); %preallocate structure (one per block)
Bmx=cell(n_mx_blocks,1); %preallocate B structure (one per block)
for b=1:n_mx_blocks
    Bmx{b}=getrow(Dmx,Dmx.BN==b+n_mx_blocks*(d-1));
    
    %add untrained sequences (10 trials) and mix trained/untrained trials within block
    iStart=(nTrialsPerBlock+1+(n_untrained_trials_per_block*(b-1))); iEnd=(nTrialsPerBlock+n_untrained_trials_per_block*b);
    Bmx{b}.seqNum(1:n_untrained_trials_per_block,1)=Dun.seqNum(iStart:iEnd);
    Bmx{b}.train(1:n_untrained_trials_per_block,1)=zeros(n_untrained_trials_per_block,1);    
    randIdx=randperm(numel(Bmx{b}.seqNum));%randomize indices for seqNum order
    Bmx{b}=getrow(Bmx{b},randIdx); %implement randomization
    
    %add exp timing info
    if fixedTrialDuration==1
        tEnd=numel(Bmx{b}.seqNum)*(trialDur+itiDur);
        Bmx{b}.tStart=(tStart:(trialDur+itiDur):tEnd-(trialDur+itiDur))'; %start time for each trial
    else
        Bmx{b}.tStart=ones(numel(Bmx{b}.seqNum),1)*(-1); %use invalid cue -1 to
    end
    Bmx{b}.iti=ones(numel(Bmx{b}.seqNum),1)*itiDur; %fixed inter-trial-interval
    Bmx{b}.feedback=ones(numel(Bmx{b}.seqNum),1); %give feedback (yes/no)
    Bmx{b}.hand=ones(numel(Bmx{b}.seqNum),1)*2; %right hand (left=1/right=2)
    
    %add finger press and cue press (cueP) info (specific to each seqNum)
    for t=1:numel(Bmx{b}.seqNum)
        if Bmx{b}.seqNum(t)<=nSeq_train 
            %trained seq trials
            Bmx{b}.press1(t,1)=train_set{Bmx{b}.seqNum(t)}(1);
            Bmx{b}.press2(t,1)=train_set{Bmx{b}.seqNum(t)}(2);
            Bmx{b}.press3(t,1)=train_set{Bmx{b}.seqNum(t)}(3);
            Bmx{b}.press4(t,1)=train_set{Bmx{b}.seqNum(t)}(4);
            Bmx{b}.press5(t,1)=train_set{Bmx{b}.seqNum(t)}(5);
            Bmx{b}.cueP{t,1}=num2str(train_set{Bmx{b}.seqNum(t)},'%d%d%d%d%d');
        else
            %untrained seq trials
            Bmx{b}.press1(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,1);
            Bmx{b}.press2(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,2);
            Bmx{b}.press3(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,3);
            Bmx{b}.press4(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,4);
            Bmx{b}.press5(t,1)=untrained_set(Bmx{b}.seqNum(t)-nSeq_train,5);
            Bmx{b}.cueP{t,1}=num2str(untrained_set(Bmx{b}.seqNum(t)-nSeq_train,:),'%d%d%d%d%d');
        end
    end
    
    %% remove extra fields (not needed for actual target files)
    extraFields={'BN'};
    Bmx{b}=rmfield(Bmx{b},extraFields);
    
    %% save structure Bmx{b} as a target file (.tgt)
    filename.fn{b,1}=fullfile(targetFolder,sprintf('sr2_fc_mx_s%02d_b%02d.tgt',s,b+n_mx_blocks*(d-1)));
    dsave(filename.fn{b,1},Bmx{b})
end

%% return output data structure D
varargout{1}=filename;
varargout{2}=Btr;
varargout{3}=Bun;
varargout{4}=Bmx;
end