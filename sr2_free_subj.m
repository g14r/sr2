function S=sr2_free_subj(subj,fig,block,trial)
%% function S=sr2_free_subj(subj,fig,block,trial)
% Subject routine for SequenceRepetition experiment (sr2_free), free variant
%
% Example calls:
%               S=sr2_free_subj('s01');         % which subj 
%               S=sr2_free_subj('s01',0);       % which subj, with/without plot
%               S=sr2_free_subj('s01',0,3);     % which subj, with/without plot, which block, all trials
%               S=sr2_free_subj('s04',1,8,3);   % which subj, with/without plot, which block, and which trial in this block
%
%%
if nargin<2
    fig=0; %don't produce a plot
end;
%%
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr2'; %path to data
datafilename=fullfile(pathToData,sprintf('sr2_%s_free.dat',subj)); %input
outfilename=fullfile(pathToData,sprintf('sr2_%s_free.mat',subj)); %output
%%
S=[]; %preallocate an empty output structure
D=dload(datafilename); %load dataset for this subj
if (nargin<3)
    trials=1:numel(D.TN); %analyze all trials for all blocks
    block=-1; %initialize block count variable
else
    if (nargin<4)
        trials=find(D.BN==block & D.TN==1):find(D.BN==block & D.TN==numel(unique(D.TN))); %analyze all trials for this block
    else
        trials=find(D.BN==block & D.TN==trial); %analyze this trial of this block
    end;
end;
%%
for t=trials
    if ~exist('MOV','var') || (D.BN(t)~=block)
        block=D.BN(t); %update block number within the loop
        MOV=movload(fullfile(pathToData,sprintf('sr2_%s_free_%02d.mov',subj,block))); %load MOV file for this block
    end;
    fprintf(1,'\nsubject: %s_free   block: %02d   trial: %02d\n\n',subj,block,D.TN(t));
    fig_name=sprintf('sr2_free_%s_b%02d_t%02d',subj,block,D.TN(t));
    T=sr2_free_trial(MOV{D.TN(t)},getrow(D,t),fig,fig_name);
    S=addstruct(S,T,'row','force');
end
%%
if (nargin<3)
    save(outfilename,'-struct','S');
end;
