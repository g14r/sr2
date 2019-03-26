%% sr2_training

% paths
pathToAnalyze='/Users/gariani/Documents/data/SequenceRepetition/sr2/analyze';

% subjects
subj={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20'};
ns=numel(subj);
subvec=zeros(1,ns);
for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end

% load group data
D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));

% add IPI info
D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);

% IPIs Comparison
T = tapply(D,{'SN', 'train'},...
    {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
    {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
    'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10 & D.prepTime==2400);
T2.SN = [T.SN; T.SN];
T2.IPI = [T.first2ipi; T.last2ipi];
T2.IPIpair = [ones(numel(T.first2ipi),1); ones(numel(T.last2ipi),1)*2];
T2.train = [T.train; T.train];

% stats
pivottable(T2.SN, [T2.train T2.IPIpair], T2.IPI, 'numel');
T.anova = anovaMixed(T2.IPI, T2.SN, 'within', [T2.IPIpair, T2.train], {'IPIpair','train'});





%% sr2_free

% paths
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr2';
pathToAnalyze='/Users/gariani/Documents/data/SequenceRepetition/sr2/analyze';
pathToDropbox='/Users/gariani/Dropbox/sr2';

% subjects
subj={'s01','s02','s03','s04','s05','s06','s08','s09','s11','s13','s15','s16','s18','s19','s20'};
ns=numel(subj);
subvec=zeros(1,ns);
for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end

% load group data
D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));

% add IPI info
D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);

% IPIs Compariso
T = tapply(D,{'SN', 'train'},...
    {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
    {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
    'subset',D.isError==0 & D.free==0 & D.prepTime==2400);
T3.SN = [T.SN; T.SN];
T3.IPI = [T.first2ipi; T.last2ipi];
T3.IPIpair = [ones(numel(T.first2ipi),1); ones(numel(T.last2ipi),1)*2];
T3.train = [T.train; T.train];

% stats
pivottable(T3.SN, [T3.train T3.IPIpair], T3.IPI, 'numel');
T.anova = anovaMixed(T3.IPI, T3.SN, 'within', [T3.IPIpair, T3.train], {'IPIpair','train'});




T4 = addstruct(T2, T3);

% stats
pivottable(T4.SN, [T4.train T4.IPIpair], T4.IPI, 'numel');
T.anova = anovaMixed(T4.IPI, T4.SN, 'within', [T4.IPIpair, T4.train], {'IPIpair','train'});
