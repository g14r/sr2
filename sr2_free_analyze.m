function[varargout]=sr2_free_analyze(what,varargin)
%% [varargout]=sr2_free_analyze(what,varargin)
% SequenceRepetition experiment (sr2_free) analysis of behavioral data
%
% usage|example calls:
%
%                       sr2_free_analyze('all_subj');                               %pre-analysis: run the subject routine for all_subj
%                       sr2_free_analyze('all_subj',{'s01'});                       %pre-analysis: run the subject routine for selected subjects
%                       [all_data]=sr2_free_analyze('all_data');                    %pre-analysis: create .mat file containing data from subjects in subj
%
%                       [D]=sr2_free_analyze('task');                               %group          analysis of task (free/forced RT) effects over time on MT and ER
%                       [D]=sr2_free_analyze('task','s01');                         %single subject analysis of task (free/forced RT) effects over time on MT and ER
%
%                       [D]=sr2_free_analyze('prepTime');                           %group          analysis of the effects of prep time on MT and ER, for trained/random sequences
%                       [D]=sr2_free_analyze('prepTime','s01');                     %single subject analysis of the effects of prep time on MT and ER, for trained/random sequences
%
%                       [D]=sr2_free_analyze('RT');                                 %group          analysis of the distribution of reaction times
%                       [D]=sr2_free_analyze('RT','s01');                           %single subject analysis of the distribution of reaction times
%
%                       [D]=sr2_free_analyze('IPI');                                %group          analysis of the inter-press intervals (IPIs) as a function of training and prep time
%                       [D]=sr2_free_analyze('IPI','s01');                          %single subject analysis of the inter-press intervals (IPIs) as a function of training and prep time
%
%                       [D]=sr2_free_analyze('memory');                             %group          analysis of memory tests: free sequence recall and sequence recognition (pre/post sequence production)
%                       [D]=sr2_free_analyze('memory','s01');                       %single subject analysis of memory tests: free sequence recall and sequence recognition (pre/post sequence production)
%
%                       [D]=sr2_free_analyze('date');                               %group          analysis of time difference between day 4 (test) and follow-up session (retest)
%                       [D]=sr2_free_analyze('date','s01');                         %single subject analysis of time difference between day 4 (test) and follow-up session (retest)
%
%                       [D]=sr2_free_analyze('force');                              %group          data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
%                       [D]=sr2_free_analyze('force','s01');                        %single subject data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
%
%                       [D]=sr2_free_analyze('peak_force');                         %group          peak press force
%                       [D]=sr2_free_analyze('peak_force','s01');                   %single subject peak press force
%
%                       [D]=sr2_free_analyze('press_dur_ms');                       %group          analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_free_analyze('press_dur_ms','s01');                 %single subject analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%
%
% --
% gariani@uwo.ca - 2018.01.24

%% globals

% paths
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr2';
pathToAnalyze='/Users/gariani/Documents/data/SequenceRepetition/sr2/analyze';
pathToDropbox='/Users/gariani/Dropbox/sr2';

% subjects
subj={'s01','s02','s03','s04','s05','s06','s08','s09','s11','s13','s15','s16','s18','s19','s20'};
ns=numel(subj);
subvec=zeros(1,ns);
for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue=[49,130,189]/255;
lightblue=[158,202,225]/255;
red=[222,45,38]/255;
lightred=[252,146,114]/255;
green=[49,163,84]/255;
lightgreen=[161,217,155]/255;
orange=[253,141,60]/255;
yellow=[254,196,79]/255;
lightyellow=[255,237,160]/255;
purple=[117,107,177]/255;
lightpurple=[188,189,220]/255;
darkgray=[50,50,50]/255;
gray=[150,150,150]/255;
lightgray=[189,189,189]/255;
silver=[240,240,240]/255;
black=[0,0,0]/255;
ptc_un_cbs={[229 158 102]/255;  [221 126 50]/255;  [213 94 0]/255;   [149 65 0]/255;  [50 50 50]/255};
ptc_tr_cbs={[153 198 224]/255;  [76 156 201]/255;  [0 102 160]/255;  [0 68 106]/255;  [50 50 50]/255};

% plot defaults
fs=20; %default font size for all figures
lw=4; %3; %default line width for all figures
ms=13; %10; %default marker size for all figures

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
%trsty=style.custom({red,green},'linewidth',lw,'markersize',ms,'errorwidth',lw);
%trsty_cbs=style.custom({cbs_red,cbs_blue},'linewidth',lw,'markersize',ms,'errorwidth',lw);
ptUNsty=style.custom(ptc_un_cbs,'markersize',ms,'linewidth',lw);
ptTRsty=style.custom(ptc_tr_cbs,'markersize',ms,'linewidth',lw);
tfsty_cbs=style.custom({cbs_red,cbs_blue,gray,darkgray},'linewidth',lw,'markersize',ms,'errorwidth',lw);
dursty=style.custom({lightblue,blue},'markersize',ms,'linewidth',lw);
delsty=style.custom({lightred,red},'markersize',ms,'linewidth',lw);

% legends
sessleg={'Pre SP','Post SP','Avg both'};
trleg={'Untrained','Trained'};
ptleg={'400','800','1600','2400','free-RT'};
tfleg={'Untrained','Trained','Untrained free-RT','Trained free-RT'};

% hands
LH=[1 2 3 4 5]; %left hand column indices (force analysis)
RH=[6 7 8 9 10]; %right hand column indices (force analysis)

%% types of analysis
switch (what)
    case 'all_subj' % pre-analysis: run the subject routine for all_subj
        if nargin>1; subj=varargin{1}; end
        for s=1:numel(subj)
            sr2_free_subj(subj{s},0); % run sr2_free_subj.m routine (without plot)
        end
        
    case 'all_data' % pre-analysis: create .mat file containing data from all subjects
        all_data=[];
        for s=1:ns
            fprintf('\n%s\n\n',subj{s});
            D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj{s}))); %load data structure for each subject
            D.SN=ones(numel(D.TN),1)*subvec(s); % add information about subject number
            all_data=addstruct(all_data,D); % append data structures from each subject
        end
        save(fullfile(pathToAnalyze,'sr2_free_all_data.mat'),'-struct', 'all_data'); %save all_data.mat file
        varargout={all_data};
        
    case 'task' % analysis of training effects over time on MT and ER, for different prep times
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        end
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        D.firstIPIs=sum([D.IPI_1,D.IPI_2],2); D.lastIPIs=sum([D.IPI_3,D.IPI_4],2);
        
        % create summary tables for MT, timing ER, press ER
        D1=tapply(D,{'SN','day','train'},...
            {D.MT,'nanmean','name','MTfc','subset',D.free==0},...
            {D.MT,'nanmean','name','MTfr','subset',D.free==1},...
            {D.MT,'nanmean','name','MTb','subset',D.mix==0},...
            {D.MT,'nanmean','name','MTm','subset',D.mix==1},...
            {D.MT,'nanmean','name','MTfcb','subset',D.free==0 & D.mix==0},...
            {D.MT,'nanmean','name','MTfrb','subset',D.free==1 & D.mix==0},...
            {D.MT,'nanmean','name','MTfcm','subset',D.free==0 & D.mix==1},...
            {D.MT,'nanmean','name','MTfrm','subset',D.free==1 & D.mix==1},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','day','train'},...
            {D.timingError,'nanmean','name','timingErrorfc','subset',D.free==0},...
            {D.timingError,'nanmean','name','timingErrorfr','subset',D.free==1},...
            {D.timingError,'nanmean','name','timingErrorb','subset',D.mix==0},...
            {D.timingError,'nanmean','name','timingErrorm','subset',D.mix==1},...
            {D.timingError,'nanmean','name','timingErrorfcb','subset',D.free==0 & D.mix==0},...
            {D.timingError,'nanmean','name','timingErrorfrb','subset',D.free==1 & D.mix==0},...
            {D.timingError,'nanmean','name','timingErrorfcm','subset',D.free==0 & D.mix==1},...
            {D.timingError,'nanmean','name','timingErrorfrm','subset',D.free==1 & D.mix==1},...
            'subset',D.dummy==0 & D.rtt==0);
        D3=tapply(D,{'SN','day','train'},...
            {D.pressError,'nanmean','name','pressErrorfc','subset',D.free==0},...
            {D.pressError,'nanmean','name','pressErrorfr','subset',D.free==1},...
            {D.pressError,'nanmean','name','pressErrorb','subset',D.mix==0},...
            {D.pressError,'nanmean','name','pressErrorm','subset',D.mix==1},...
            {D.pressError,'nanmean','name','pressErrorfcb','subset',D.free==0 & D.mix==0},...
            {D.pressError,'nanmean','name','pressErrorfrb','subset',D.free==1 & D.mix==0},...
            {D.pressError,'nanmean','name','pressErrorfcm','subset',D.free==0 & D.mix==1},...
            {D.pressError,'nanmean','name','pressErrorfrm','subset',D.free==1 & D.mix==1},...
            'subset',D.timingError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Task effects - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Task effects - group (N=%d)',ns)); end
        
        % MT
        sty=style.custom({lightblue,blue},'markersize',ms,'linewidth',lw);
        subplot(3,8,1);
        plt.bar(D1.day,D1.MTfc,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Forced-RT');
        subplot(3,8,2);
        plt.bar(D1.day,D1.MTfcb,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Forced-RT | Blocked');
        subplot(3,8,3);
        plt.bar(D1.day,D1.MTfcm,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Forced-RT | Mixed');
        subplot(3,8,4);
        plt.bar(D1.day,D1.MTfr,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Free-RT');
        subplot(3,8,5);
        plt.bar(D1.day,D1.MTfrb,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Free-RT | Blocked');
        subplot(3,8,6);
        plt.bar(D1.day,D1.MTfrm,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Free-RT | Mixed');
        subplot(3,8,7);
        plt.bar(D1.day,D1.MTb,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Blocked');
        subplot(3,8,8);
        plt.bar(D1.day,D1.MTm,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); ylim([700 1700]); title('Mixed');
        
        % tER
        sty=style.custom({lightred,red},'markersize',ms,'linewidth',lw);
        subplot(3,8,9);
        plt.bar(D2.day,D2.timingErrorfc*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,10);
        plt.bar(D2.day,D2.timingErrorfcb*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,11);
        plt.bar(D2.day,D2.timingErrorfcm*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,12);
        plt.bar(D2.day,D2.timingErrorfr*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,13);
        plt.bar(D2.day,D2.timingErrorfrb*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,14);
        plt.bar(D2.day,D2.timingErrorfrm*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,15);
        plt.bar(D2.day,D2.timingErrorb*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        subplot(3,8,16);
        plt.bar(D2.day,D2.timingErrorm*100,'split',D2.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Timing error (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100);
        
        % pER
        sty=style.custom({lightgreen,green},'markersize',ms,'linewidth',lw);
        subplot(3,8,17);
        plt.bar(D3.day,(1-D3.pressErrorfc)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,18);
        plt.bar(D3.day,(1-D3.pressErrorfcb)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,19);
        plt.bar(D3.day,(1-D3.pressErrorfcm)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,20);
        plt.bar(D3.day,(1-D3.pressErrorfr)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,21);
        plt.bar(D3.day,(1-D3.pressErrorfrb)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,22);
        plt.bar(D3.day,(1-D3.pressErrorfrm)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,23);
        plt.bar(D3.day,(1-D3.pressErrorb)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        subplot(3,8,24);
        plt.bar(D3.day,(1-D3.pressErrorm)*100,'split',D3.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); yticks(0:20:100); drawline(20,'dir','horz','linestyle','--');
        
        % out
        D.D1=D1; D.D2=D2; D.D3=D3; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'prepTime' % analysis of the effects of prep time on MT and ER, for trained/random sequences
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        end
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        
        % compute actual preparation time for forced-RT trials
        D.realPrepTime(D.free==0,1)=D.prepTime(D.free==0,1)+D.RT(D.free==0,1);
        D.realPrepTime(D.free==1,1)=D.RT(D.free==1,1);
        
        % divide free-RT in quartiles and assign them to prepTime categories
        D = getrow(D, D.isError==0);
        [i,~,~] = split_data(D.realPrepTime, 'split',[D.SN D.train D.mix], 'numquant',3, 'subset',D.free==1);
        D.prepTime(D.free==1 & i==1) = 400;
        D.prepTime(D.free==1 & i==2) = 800;
        D.prepTime(D.free==1 & i==3) = 1600;
        D.prepTime(D.free==1 & i==4) = 2400;
        
        % create summary table for MT
        T = tapply(D, {'SN','train','prepTime','mix','free'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.realPrepTime,'nanmean', 'name','realPrepTime'}, ...
            'subset',D.isError==0);
        
%         % divide free-RT in quartiles and assign them to prepTime categories
%         [i,~,~] = split_data(T.realPrepTime, 'split',[T.train T.mix], 'numquant',4, 'subset',T.free==1);
%         T.prepTime(T.free==1 & i==1) = 400;
%         T.prepTime(T.free==1 & i==2) = 800;
%         T.prepTime(T.free==1 & i==3) = 1600;
%         T.prepTime(T.free==1 & i==4) = 2400;

        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        %         % compute the percentage of MT at prep time 2400 to show the interaction
        %         h=figure; [~,y]=xyplot(T.realPrepTime,(1./T.normMT),T.prepTime,'split',[T.train T.free T.mix],'style_thickline','leg','auto'); close(h);
        %         T.pctMT = zeros(numel(T.MT),1);
        %         trbidx = T.train==1 & T.mix==0 & T.free==0;
        %         unbidx = T.train==0 & T.mix==0 & T.free==0;
        %         trmidx = T.train==1 & T.mix==1 & T.free==0;
        %         unmidx = T.train==0 & T.mix==1 & T.free==0;
        %         %         T.pctMT(unbidx) = 100 + (100 - ( T.normMT(unbidx) * 100 ) / y(1,end));
        %         %         T.pctMT(unmidx) = 100 + (100 - ( T.normMT(unmidx) * 100 ) / y(2,end));
        %         %         T.pctMT(trbidx) = 100 + (100 - ( T.normMT(trbidx) * 100 ) / y(5,end));
        %         %         T.pctMT(trmidx) = 100 + (100 - ( T.normMT(trmidx) * 100 ) / y(6,end));
        %         T.pctMT(unbidx) = (1 ./ (T.normMT(unbidx)))  ./  y(1,end) * 100;
        %         T.pctMT(unmidx) = (1 ./ (T.normMT(unmidx)))  ./  y(2,end) * 100;
        %         T.pctMT(trbidx) = (1 ./ (T.normMT(trbidx)))  ./  y(5,end) * 100;
        %         T.pctMT(trmidx) = (1 ./ (T.normMT(trmidx)))  ./  y(6,end) * 100;
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.train], [T.prepTime,T.mix], T.MT, 'length');

        % open figure
        if nargin>1; figure('Name',sprintf('Prep time - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Prep time - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % plot data
        subplot(2,2,1); title('Blocked');
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'split',[T.free T.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','north', 'subset',T.mix==0);
        xlabel('Preparation time (ms)'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(T.free==0))); 
        xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([840 1610]);
        hold on; drawline(unique(T.prepTime(T.free==0)),'dir','vert','linestyle',':','color','k'); hold off;
        
        %         subplot(2,2,2); title('Blocked, Interaction');
        %         plt.xy(T.realPrepTime,T.pctMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','east', 'subset',T.free==0 & T.mix==0);
        %         xlabel('Preparation time (ms)'); ylabel('% of ET at 2400'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(T.free==0))); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([65 105]);
        %         hold on; drawline(unique(T.prepTime(T.free==0)),'dir','vert','linestyle',':','color','k'); drawline(100,'dir','horz','linestyle','--','color','k'); hold off;
        %
        subplot(2,2,2); title('Mixed');
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'split',[T.free T.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','north', 'subset',T.mix==1);
        xlabel('Preparation time (ms)'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(T.free==0))); 
        xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([840 1610]);
        hold on; drawline(unique(T.prepTime(T.free==0)),'dir','vert','linestyle',':','color','k'); hold off;
        
        %         subplot(2,2,4); title('Mixed, Interaction');
        %         plt.xy(T.realPrepTime,T.pctMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','east', 'subset',T.free==0 & T.mix==1);
        %         xlabel('Preparation time (ms)'); ylabel('% of ET at 2400'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(T.free==0))); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([65 105]);
        %         hold on; drawline(unique(T.prepTime(T.free==0)),'dir','vert','linestyle',':','color','k'); drawline(100,'dir','horz','linestyle','--','color','k'); hold off;
        %
        
        % stats
        %         T = tapply(T, {'SN','train'}, ...
        %             {T.realPrepTime,'nanmean', 'name','RT'}, ...
        %             {T.realPrepTime,'nanmean', 'name','RTb', 'subset',T.mix==0}, ...
        %             {T.realPrepTime,'nanmean', 'name','RTm', 'subset',T.mix==1}, ...
        %             'subset',T.free==1);
        %         ttest(T.RT(T.train==0), T.RT(T.train==1), 2, 'paired');
        %         ttest(T.RTb(T.train==0), T.RTb(T.train==1), 2, 'paired');
        %         ttest(T.RTm(T.train==0), T.RTm(T.train==1), 2, 'paired');
        
        T = tapply(D, {'SN','free'}, ...
            {D.MT,'nanmean', 'name','ET', 'subset',D.prepTime==0 | D.prepTime>800}, ...
            'subset',D.isError==0);
        ttest(T.ET(T.free==0), T.ET(T.free==1), 2, 'paired');
        
        %         T = tapply(D, {'SN','train','mix'}, ...
        %             {D.MT,'nanmean', 'name','MT'}, ...
        %             'subset',D.isError==0 & D.free==1);
        %         ttest(T.MT(T.mix==0 & T.train==0), T.MT(T.mix==0 & T.train==1), 2, 'paired');
        %         ttest(T.MT(T.mix==1 & T.train==0), T.MT(T.mix==1 & T.train==1), 2, 'paired');
        
        T = tapply(D, {'SN','train'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.isError==0 & D.free==0 & D.mix==0 & D.prepTime==2400);
        ttest(T.MT(T.train==0), T.MT(T.train==1), 2, 'paired');
        
        T = tapply(D, {'SN','train'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.isError==0 & D.free==0 & D.mix==1 & D.prepTime==2400);
        ttest(T.MT(T.train==0), T.MT(T.train==1), 2, 'paired');
        
        D = getrow(D, D.prepTime>0);
        T = tapply(D, {'SN','train','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.isError==0);
        T.anova = anovaMixed(T.MT,T.SN,'within', [T.prepTime,T.train], {'prepTime','train'});
        
        %
        %         T = tapply(D, {'SN','train','prepTime'}, ...
        %             {D.MT,'nanmean', 'name','MT'}, ...
        %             'subset',D.isError==0 & D.mix==0);
        %         T.anova_blocked = anovaMixed(T.MT,T.SN,'within', [T.prepTime,T.train], {'prepTime','train'});
        %
        %         T = tapply(D, {'SN','train','prepTime'}, ...
        %             {D.MT,'nanmean', 'name','MT'}, ...
        %             'subset',D.isError==0 & D.mix==1);
        %         T.anova_mixed = anovaMixed(T.MT,T.SN,'within', [T.prepTime,T.train], {'prepTime','train'});
        %
        %         T = tapply(D, {'SN','train'}, {D.MT,'nanmean', 'name','MT'}, 'subset',D.isError==0 & D.mix==1 & D.prepTime==2400);
        %         ttest(T.MT(T.train==0), T.MT(T.train==1), 2, 'paired');
        
        %         if nargin>1 % load single subj data
        %             subj=varargin{1};
        %             D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj))); %load data for this subject
        %             D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        %         else % load group data
        %             D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        %         end
        %
        %         % add IPI info
        %         D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        %         D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        %
        %         % compute actual preparation time for forced-RT trials
        %         D.realPrepTime(D.free==0,1)=D.prepTime(D.free==0,1)+D.RT(D.free==0,1);
        %         D.realPrepTime(D.free==1,1)=D.RT(D.free==1,1);
        %
        %         % create summary table for MT
        %         T = tapply(D, {'SN','train','prepTime','mix','free'}, ...
        %             {D.realPrepTime,'nanmean', 'name','RPT'}, ...
        %             'subset',D.isError==0);
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T = normData(T, {'RPT'}, 'sub');
        %
        %         % open figure
        %         if nargin>1; figure('Name',sprintf('Training - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Training - group (N=%d)',ns)); end
        %         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        %
        %         subplot(2,2,1); title('Blocked');
        %         plt.box(T.prepTime,T.RPT, 'split',[T.free T.train], 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','north', 'subset',T.mix==0);
        %         xlabel('Preparation time (ms)'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        %         xticks(unique(T.prepTime)); %xlim([min(T.prepTime)-100 max(T.prepTime)+100]); %ylim([850 1600]);
        %         hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        %
        %         subplot(2,2,2); title('Mixed');
        %         plt.box(T.prepTime,T.RPT, 'split',[T.free T.train], 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','north', 'subset',T.mix==1);
        %         xlabel('Preparation time (ms)'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        %         xticks(unique(T.prepTime)); %xlim([min(T.prepTime)-100 max(T.prepTime)+100]); %ylim([850 1600]);
        %         hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'RT' % analysis of the distribution of reaction times
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        end
        
        % create summary tables for MT, timing ER, press ER
        D1=tapply(D,{'SN','day','train'},...
            {D.RT,'nanmean','name','RTfc','subset',D.free==0},...
            {D.RT,'nanmean','name','RTfr','subset',D.free==1},...
            {D.RT,'nanmean','name','RTb','subset',D.mix==0},...
            {D.RT,'nanmean','name','RTm','subset',D.mix==1},...
            {D.RT,'nanmean','name','RTfcb','subset',D.free==0 & D.mix==0},...
            {D.RT,'nanmean','name','RTfrb','subset',D.free==1 & D.mix==0},...
            {D.RT,'nanmean','name','RTfcm','subset',D.free==0 & D.mix==1},...
            {D.RT,'nanmean','name','RTfrm','subset',D.free==1 & D.mix==1},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        x=500; %define x-axis range/cut-off
        D2=getrow(D,D.dummy==0 & D.rtt==0 & D.RT<=x & D.RT>=-x);
        
        D3=tapply(D,{'SN','BN','train'},...
            {D.RT,'nanmean','name','RTfc','subset',D.free==0},...
            {D.RT,'nanmean','name','RTfr','subset',D.free==1},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        if nargin>1; figure('Name',sprintf('RT analysis - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('RT analysis - group (N=%d)',ns)); end
        sty=style.custom({lightblue,blue},'markersize',ms,'linewidth',lw);
        
        subplot(3,6,1);
        nurt=numel(unique(D2.RT(D2.free==0)));
        title('Forced-RT'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]); axis square;
        plt.hist(D2.RT,'split',D2.train,'numcat',nurt,'percent',1,'style',sty,'subset',D2.free==0,'leg',trleg); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        
        subplot(3,6,2);
        nurt=numel(unique(D2.RT(D2.free==0 & D2.mix==0)));
        title('Forced-RT | Blocked'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]); axis square;
        plt.hist(D2.RT,'split',D2.train,'numcat',nurt,'percent',1,'style',sty,'subset',D2.free==0 & D2.mix==0,'leg',trleg); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        
        subplot(3,6,3);
        nurt=numel(unique(D2.RT(D2.free==0 & D2.mix==1)));
        title('Forced-RT | Mixed'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]); axis square;
        plt.hist(D2.RT,'split',D2.train,'numcat',nurt,'percent',1,'style',sty,'subset',D2.free==0 & D2.mix==1,'leg',trleg); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        
        subplot(3,6,4);
        nurt=numel(unique(D.RT(D.free==1)));
        title('Free-RT'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); axis square;
        plt.hist(D.RT,'split',D.train,'numcat',nurt,'percent',1,'style',sty,'subset',D.free==1,'leg',trleg);
        
        subplot(3,6,5);
        nurt=numel(unique(D.RT(D.free==1 & D.mix==0)));
        title('Free-RT | Blocked'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); axis square;
        plt.hist(D.RT,'split',D.train,'numcat',nurt,'percent',1,'style',sty,'subset',D.free==1 & D.mix==0,'leg',trleg);
        
        subplot(3,6,6);
        nurt=numel(unique(D.RT(D.free==1 & D.mix==1)));
        title('Free-RT | Mixed'); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); axis square;
        plt.hist(D.RT,'split',D.train,'numcat',nurt,'percent',1,'style',sty,'subset',D.free==1 & D.mix==1,'leg',trleg);
        
        subplot(3,6,7);
        plt.bar(D1.day,D1.RTfc,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        subplot(3,6,8);
        plt.bar(D1.day,D1.RTfcb,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        subplot(3,6,9);
        plt.bar(D1.day,D1.RTfcm,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square;
        
        subplot(3,6,10);
        plt.bar(D1.day,D1.RTfr,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([1000 1600]);
        
        subplot(3,6,11);
        plt.bar(D1.day,D1.RTfrb,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([1000 1600]);
        
        subplot(3,6,12);
        plt.bar(D1.day,D1.RTfrm,'split',D1.train,'style',sty,'leg',trleg);
        xticklabels(''); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); axis square; ylim([1000 1600]);
        
        subplot(3,6,13:15);
        plt.bar(D3.SN(~isnan(D3.RTfc)),D3.RTfc(~isnan(D3.RTfc)),'split',D3.train,'style',sty,'leg',trleg);
        xlabel('Subject number'); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); title('Forced-RT');
        
        subplot(3,6,16:18);
        plt.bar(D3.SN(~isnan(D3.RTfr)),D3.RTfr(~isnan(D3.RTfr)),'split',D3.train,'style',sty,'leg',trleg);
        xlabel('Subject number'); ylabel('Mean RT (ms)'); set(gca,'fontsize',fs); title('Free-RT');
        
        %         % Scatter plot MT-RT
        %         if nargin>1; figure('Name',sprintf('Scatterplot MT-RT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Scatterplot RT-MT - group (N=%d)',ns)); end
        %         sty=style.custom({red,blue},'markersize',ms,'linewidth',lw);
        %
        %         subplot(1,3,1);
        %         plt.scatter(D.MT,D.RT,'split',D.train,'subset',D.free==1 & D.dummy==0 & D.isError==0 & D.rtt==0,'style',sty,'leg',trleg);
        %         title('Free-RT'); xlabel('MT (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs); %axis square;
        %
        %         subplot(1,3,2);
        %         plt.scatter(D.MT,D.RT,'split',D.train,'subset',D.free==1 & D.mix==0 & D.dummy==0 & D.isError==0 & D.rtt==0,'style',sty,'leg',trleg);
        %         title('Free-RT | Blocked'); xlabel('MT (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs); %axis square;
        %
        %         subplot(1,3,3);
        %         plt.scatter(D.MT,D.RT,'split',D.train,'subset',D.free==1 & D.mix==1 & D.dummy==0 & D.isError==0 & D.rtt==0,'style',sty,'leg',trleg);
        %         title('Free-RT | Mixed'); xlabel('MT (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs); %axis square;
        
        % stats
        T = tapply(D,{'SN','train'},...
            {D.RT,'nanmean','name','RTfr','subset',D.free==1},...
            {D.RT,'nanmean','name','RTfrb','subset',D.free==1 & D.mix==0},...
            {D.RT,'nanmean','name','RTfrm','subset',D.free==1 & D.mix==1},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        ttest(T.RTfr(T.train==0), T.RTfr(T.train==1), 2, 'paired');
        ttest(T.RTfrb(T.train==0), T.RTfrb(T.train==1), 2, 'paired');
        ttest(T.RTfrm(T.train==0), T.RTfrm(T.train==1), 2, 'paired');
        
        % out
        D.D1=D1; D.D2=D2; D.D3=D3; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'IPI' % analysis of the inter-press intervals (IPIs) as a function of training and prep time
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s_free.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        end
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        
        % create summary table for IPIs
        T=tapply(D,{'SN','prepTime','train','mix','free'},...
            {D.IPI_1,'nanmean','name','IPI1'},...
            {D.IPI_2,'nanmean','name','IPI2'},...
            {D.IPI_3,'nanmean','name','IPI3'},...
            {D.IPI_4,'nanmean','name','IPI4'},...
            'subset',D.isError==0);
        for i=1:size(D.IPI,2)
            T.IPI(:,i)=eval(sprintf('T.IPI%d',i));
            T=rmfield(T,sprintf('IPI%d',i));
            T.IPInum(:,i)=ones(size(T.SN,1),1)*i;
            T.SN(:,i)=T.SN(:,1);
            T.prepTime(:,i)=T.prepTime(:,1);
            T.train(:,i)=T.train(:,1);
            T.mix(:,i)=T.mix(:,1);
            T.free(:,i)=T.free(:,1);
        end
        T.IPI=reshape(T.IPI,size(T.IPI,1)*size(T.IPI,2),1);
        T.IPInum=reshape(T.IPInum,size(T.IPInum,1)*size(T.IPInum,2),1);
        T.SN=reshape(T.SN,size(T.SN,1)*size(T.SN,2),1);
        T.prepTime=reshape(T.prepTime,size(T.prepTime,1)*size(T.prepTime,2),1);
        T.train=reshape(T.train,size(T.train,1)*size(T.train,2),1);
        T.mix=reshape(T.mix,size(T.mix,1)*size(T.mix,2),1);
        T.free=reshape(T.free,size(T.free,1)*size(T.free,2),1);
        
        % normalize IPI data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable([T.train T.IPInum], [T.mix T.prepTime], T.normIPI, 'length');
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1); title('Untrained, Blocked');
        plt.line(T.IPInum,T.normIPI, 'split',[T.free T.prepTime], 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==0 & T.mix==0);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        subplot(2,2,2); title('Trained, Blocked');
        plt.line(T.IPInum,T.normIPI, 'split',[T.free T.prepTime], 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==1 & T.mix==0);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        subplot(2,2,3); title('Untrained, Mixed');
        plt.line(T.IPInum,T.normIPI, 'split',[T.free T.prepTime], 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==0 & T.mix==1);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        subplot(2,2,4); title('Trained, Mixed');
        plt.line(T.IPInum,T.normIPI, 'split',[T.free T.prepTime], 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==1 & T.mix==1);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        %         % stats
        %         fprintf(1, '\nANOVA:\n');
        %         anovaMixed(T.IPI, T.SN, 'within',[T.prepTime, T.train, T.mix],{'prepTime','train','mix'}, 'subset',T.free==0 & T.IPInum==4);
        %
        %         % trained blocked
        %         fprintf(1, '\ntrained blocked:\n');
        %         ttest(T.IPI(T.prepTime==400 & ismember(T.IPInum,4) & T.train==1 & T.free==0 & T.mix==0), T.IPI(T.prepTime==2400 & ismember(T.IPInum,4) & T.train==1 & T.free==0 & T.mix==0), 2, 'paired');
        %
        %         % trained mixed
        %         fprintf(1, '\ntrained mixed:\n');
        %         ttest(T.IPI(T.prepTime==400 & ismember(T.IPInum,4) & T.train==1 & T.free==0 & T.mix==1), T.IPI(T.prepTime==2400 & ismember(T.IPInum,4) & T.train==1 & T.free==0 & T.mix==1), 2, 'paired');
        %
        %         % untrained blocked
        %         fprintf(1, '\nuntrained blocked:\n');
        %         ttest(T.IPI(T.prepTime==400 & ismember(T.IPInum,4) & T.train==0 & T.free==0 & T.mix==0), T.IPI(T.prepTime==2400 & ismember(T.IPInum,4) & T.train==0 & T.free==0 & T.mix==0), 2, 'paired');
        %
        %         % untrained mixed
        %         fprintf(1, '\nuntrained mixed:\n');
        %         ttest(T.IPI(T.prepTime==400 & ismember(T.IPInum,4) & T.train==0 & T.free==0 & T.mix==1), T.IPI(T.prepTime==2400 & ismember(T.IPInum,4) & T.train==0 & T.free==0 & T.mix==1), 2, 'paired');
        
        %---------------------------------------------------------------------------------------------------
        % IPIs Comparison
        % compute actual preparation time for forced-RT trials
        D.realPrepTime(D.free==0,1)=D.prepTime(D.free==0,1)+D.RT(D.free==0,1);
        D.realPrepTime(D.free==1,1)=D.RT(D.free==1,1);
        
        % divide free-RT in quartiles and assign them to prepTime categories
        [i,~,~] = split_data(D.realPrepTime, 'split',[D.SN D.train D.mix], 'numquant',4, 'subset',D.isError==0 & D.free==1);
        D.prepTime(D.free==1 & i==1) = 400;
        D.prepTime(D.free==1 & i==2) = 800;
        D.prepTime(D.free==1 & i==3) = 1600;
        D.prepTime(D.free==1 & i==4) = 2400;
        
        % create summary table
        T2=tapply(D,{'SN','prepTime','train','mix','free'},...
            {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...            %{D.IPI_3,'nanmean','name','IPI3'},...             %{D.IPI_4,'nanmean','name','IPI4'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
            {D.realPrepTime,'nanmean', 'name','realPrepTime'}, ...
            'subset',D.isError==0);
                        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'first2ipi','last2ipi','realPrepTime'}, 'sub');
        %T2 = normData(T2, {'first2ipi', 'IPI3', 'IPI4'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T2.train, [T2.mix T2.prepTime], T2.normfirst2ipi, 'length'); % pivottable(T2.train, [T2.mix T2.prepTime], T2.normlast2ipi, 'length');
        
        gap = 2500; % how much space between IPI pair plots (should be > 2400)
        T3.IPI = [T2.normfirst2ipi; T2.normlast2ipi];
        T3.RPT = [T2.normrealPrepTime; T2.normrealPrepTime+gap];
        T3.PT = [T2.prepTime; T2.prepTime+gap];
        T3.mix = [T2.mix; T2.mix];
        T3.free = [T2.free; T2.free];
        T3.train = [T2.train; T2.train];
        T3.pair = [zeros(numel(T2.SN),1); ones(numel(T2.SN),1)] + 1;
%         T3.IPI = [T2.normfirst2ipi; T2.normIPI3; T2.normIPI4];
%         T3.RPT = [T2.realPrepTime; T2.realPrepTime+gap; T2.realPrepTime+gap*2];
%         T3.PT = [T2.prepTime; T2.prepTime+gap; T2.prepTime+gap*2];
%         T3.mix = [T2.mix; T2.mix; T2.mix];
%         T3.free = [T2.free; T2.free; T2.free];
%         T3.train = [T2.train; T2.train; T2.train];
%         T3.pair = [zeros(numel(T2.SN),1); ones(numel(T2.SN),1); ones(numel(T2.SN),1)*2] + 1;
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI comparison - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI comparison - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,3); title('Blocked'); %hold on;
        plt.xy(T3.RPT,T3.IPI,T3.PT, 'split',[T3.pair T3.free T3.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T3.mix==0);
        xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); %axis square;
        xlim([min(T3.PT)-100 max(T3.PT)+100]);
        ylim([251 849]); %ylim([321 849]); %ylim([63 876]);
        xticks(unique(T3.PT(T3.free==0))); xticklabels(repmat([400,800,1600,2400],1,3));
        drawline(unique(T3.PT(T3.free==0)),'dir','vert','linestyle',':','color','k');
        
        subplot(2,2,4); title('Mixed'); %hold on;
        plt.xy([T3.RPT],T3.IPI,T3.PT, 'split',[T3.pair T3.free T3.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T3.mix==1);
        xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); %axis square; %xticks(unique(T.prepTime)); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([60 840]);
        xlim([min(T3.PT)-100 max(T3.PT)+100]);
        ylim([251 849]); %ylim([321 849]); %ylim([63 876]); 
        xticks(unique(T3.PT(T3.free==0))); xticklabels(repmat([400,800,1600,2400],1,3));
        drawline(unique(T3.PT(T3.free==0)),'dir','vert','linestyle',':','color','k');
        
        %         subplot(2,2,1); title('First 2 IPIs, Blocked');
        %         plt.xy(T2.realPrepTime,T2.normfirst2ipi,T2.prepTime, 'split',[T2.free T2.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T2.mix==0);
        %         xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([60 840]);
        %
        %         subplot(2,2,2); title('Last 2 IPIs, Blocked');
        %         plt.xy(T2.realPrepTime,T2.normlast2ipi,T2.prepTime, 'split',[T2.free T2.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T2.mix==0);
        %         xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([60 840]);
        %
        %         subplot(2,2,3); title('First 2 IPIs, Mixed');
        %         plt.xy(T2.realPrepTime,T2.normfirst2ipi,T2.prepTime, 'split',[T2.free T2.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T2.mix==1);
        %         xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([60 840]);
        %
        %         subplot(2,2,4); title('Last 2 IPIs, Mixed');
        %         plt.xy(T2.realPrepTime,T2.normlast2ipi,T2.prepTime, 'split',[T2.free T2.train], 'errorbars','plusminus_wocap', 'style',tfsty_cbs, 'leg',tfleg, 'leglocation','northeast', 'subset',T2.mix==1);
        %         xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime(T.free==0))-100 max(T.prepTime(T.free==0))+100]); ylim([60 840]);
        
        %plt.match('y');
        
        %         hold on;
        %         subplot(2,2,1); title('First 2 IPIs, Blocked');
        %         drawline(unique(T2.prepTime(T2.free==0)),'dir','vert','linestyle',':','color','k');
        %
        %         subplot(2,2,2); title('Last 2 IPIs, Blocked');
        %         drawline(unique(T2.prepTime(T2.free==0)),'dir','vert','linestyle',':','color','k');
        %
        %         subplot(2,2,3); title('First 2 IPIs, Mixed');
        %         drawline(unique(T2.prepTime(T2.free==0)),'dir','vert','linestyle',':','color','k');
        %
        %         subplot(2,2,4); title('Last 2 IPIs, Mixed');
        %         drawline(unique(T2.prepTime(T2.free==0)),'dir','vert','linestyle',':','color','k');
        %         hold off;
        
        % stats
        T2 = tapply(D, {'SN', 'train'},...
            {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
            'subset',D.isError==0 & D.free==0 & D.prepTime==2400);
        T3.SN = [T2.SN; T2.SN];
        T3.IPI = [T2.first2ipi; T2.last2ipi];
        T3.IPIpair = [ones(numel(T2.first2ipi),1); ones(numel(T2.last2ipi),1)*2];
        T3.train = [T2.train; T2.train];
        %         T3.prepTime = [T2.prepTime; T2.prepTime];
        %         T3.mix = [T2.mix; T2.mix];
        %         T3.free = [T2.free; T2.free];
        %                 T3 = tapply(T3, {'SN','train','IPIpair'}, ...
        %                     {T3.IPI,'nanmean', 'name','IPI'}, ...
        %                     'subset',T3.prepTime==2400 & T3.free==0);
        T.anova = anovaMixed(T3.IPI, T3.SN, 'within', [T3.IPIpair, T3.train], {'IPIpair','train'});
        
        
        T = tapply(D, {'SN'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','IPI34fr', 'subset',D.free==1},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','IPI34fc', 'subset',D.free==0 & D.prepTime==2400},...
            'subset',D.isError==0 & D.train==1);
        ttest(T.IPI34fc, T.IPI34fr, 2, 'paired');
        
        %         T3 = tapply(D,{'SN','prepTime','train'},...
        %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %             'subset',D.isError==0);
        %         T3.anova    = anovaMixed(T3.last2ipi, T3.SN, 'within', [T3.prepTime, T3.train], {'prepTime','train'});
        %
        %         T4 = tapply(D,{'SN','prepTime','train','mix'},...
        %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %             'subset',D.isError==0 & D.free==0);
        %         T4.anova    = anovaMixed(T4.last2ipi, T4.SN, 'within', [T4.prepTime, T4.train], {'prepTime','train'}, 'subset',T4.mix==0); % blocked blocks
        %         T4.anova    = anovaMixed(T4.last2ipi, T4.SN, 'within', [T4.prepTime, T4.train], {'prepTime','train'}, 'subset',T4.mix==1); % mixed blocks
        %
        %         T5 = tapply(D,{'SN','prepTime','train'},...
        %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %             'subset',D.isError==0 & D.free==0);
        %         T5.anova    = anovaMixed(T5.last2ipi, T5.SN, 'within', [T5.prepTime, T5.train], {'prepTime','train'}); % forced-rt
        %
        %         T6 = tapply(D,{'SN','train'},...
        %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %             'subset',D.isError==0 & D.free==1);
        %         ttest(T6.last2ipi(T6.train==0), T6.last2ipi(T6.train==1), 2, 'paired'); % free-rt
        
        
        
        %         T.IPIpair = ismember(T.IPInum, [3,4]);
        %         fprintf(1, '\nANOVA:\n');
        %         anovaMixed(T.IPI, T.SN, 'within',[T.train, T.mix T.IPIpair],{'train','mix','IPIpair'}, 'subset',T.free==0);
        %         %         T3=tapply(D,{'SN','train'},...
        %         %             {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
        %         %             'subset',D.isError==0 & D.mix==0 & D.free==0 & D.prepTime==2400);
        %         %         ttest(T3.first2ipi(T3.train==0), T3.first2ipi(T3.train==1), 2, 'paired');
        %         %
        %         T4=tapply(D,{'SN','train'},...
        %             {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
        %             'subset',D.isError==0 & D.mix==1 & D.free==0 & D.prepTime==2400);
        %         ttest(T4.first2ipi(T4.train==0), T4.first2ipi(T4.train==1), 2, 'paired');
        %         %
        %         %         T5=tapply(D,{'SN','train'},...
        %         %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %         %             'subset',D.isError==0 & D.mix==0 & D.free==0 & D.prepTime==2400);
        %         %         ttest(T5.last2ipi(T5.train==0), T5.last2ipi(T5.train==1), 2, 'paired');
        %         %
        %         %         T6=tapply(D,{'SN','train'},...
        %         %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %         %             'subset',D.isError==0 & D.mix==1 & D.free==0 & D.prepTime==2400);
        %         %         ttest(T6.last2ipi(T6.train==0), T6.last2ipi(T6.train==1), 2, 'paired');
        
        % out
        D.T=T; D.T2=T2; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'memory' % analysis of memory tests: free sequence recall and sequence recognition
        % load excel file with sequence recognition test
        xls_fn=fullfile(pathToDropbox,'sr2_logfile.xlsx'); xls_sheet='dPrime';
        [num,text,~]=xlsread(xls_fn,xls_sheet);
        
        % store data in dataframe structure
        D.subj=text(2:end,1);
        D.SN=num(:,1);
        D.session=num(:,2);
        D.HIT=num(:,3);
        D.FA=num(:,4);
        D.MISS=num(:,5);
        D.CR=num(:,6);
        D.N_tot=num(:,7);
        D.N_yes=num(:,8);
        D.N_no=num(:,9);
        D.HIT_rate=num(:,11);
        D.FA_rate=num(:,13);
        D.seqRecall=num(:,19);
        chance=(1/num(1,18))*100;
        
        % select subjects for the analysis
        if nargin>1; subj=varargin(1); subvec=str2double(subj{1}(2:3)); end % single subject analysis
        D=getrow(D,logical(sum(D.SN==subvec,2)));
        
        % calculate d' and bias for each subject and session
        for ss=1:numel(D.SN)
            [D.dPrime(ss,1),D.bias(ss,1)]=dprime(D.HIT(ss,1),D.FA(ss,1),D.N_yes(ss,1));
        end
        
        % create summary table for d' and bias
        D1=tapply(D,{'SN'},...
            {D.dPrime,'nanmean','name','d'},...
            {D.dPrime,'nanmean','name','dFirst','subset',D.session==1},...
            {D.dPrime,'nanmean','name','dSecond','subset',D.session==2},...
            {D.bias,'nanmean','name','b'},...
            {D.bias,'nanmean','name','bFirst','subset',D.session==1},...
            {D.bias,'nanmean','name','bSecond','subset',D.session==2},...
            {D.seqRecall,'nanmean','name','r'},...
            {D.seqRecall,'nanmean','name','rFirst','subset',D.session==1},...
            {D.seqRecall,'nanmean','name','rSecond','subset',D.session==2});
        
        % open figure
        if nargin>1; figure('Name',sprintf('Sequence memory - subj %02d',subvec)); else; figure('Name',sprintf('Sequence memory - group (N=%d)',ns)); end
        
        sty=style.custom({lightgray,gray,black},'markersize',ms,'linewidth',lw);
        %         subplot(3,2,1);
        %         plt.bar(D1.SN,[D1.rFirst,D1.rSecond,D1.r],'style',sty,'leg',sessleg); legend off; title('Sequence recall test'); drawline(chance,'dir','horz','linestyle','--','color',silver);
        %         xlabel('Subject number'); ylabel('Accuacy (%)'); xt=xticks; xticks(xt(2:3:end)); xticklabels(unique(D.SN)); set(gca,'fontsize',fs); ylim([0 100]);
        subplot(2,2,1);
        plt.bar(ones(numel(D1.SN),1),[D1.rFirst,D1.rSecond,D1.r],'style',sty,'leg',sessleg); drawline(chance,'dir','horz','linestyle','--','color',silver);
        xlabel('Group average'); ylabel('Accuacy (%)'); xticklabels(''); axis square; set(gca,'fontsize',fs); ylim([0 100]);
        
        sty=style.custom({yellow,orange,red},'markersize',ms,'linewidth',lw);
        %         subplot(3,2,3);
        %         plt.bar(D1.SN,[D1.dFirst,D1.dSecond,D1.d],'style',sty,'leg',sessleg); legend off; title('Sequence recognition test');
        %         xlabel('Subject number'); ylabel('d-prime'); xt=xticks; xticks(xt(2:3:end)); xticklabels(unique(D.SN)); set(gca,'fontsize',fs);
        subplot(2,2,2);
        plt.bar(ones(numel(D1.SN),1),[D1.dFirst,D1.dSecond,D1.d],'style',sty,'leg',sessleg);
        xlabel('Group average'); ylabel('d-prime'); xticklabels(''); axis square; set(gca,'fontsize',fs); %ylim([0 2]);
        
        %         sty=style.custom({lightblue,blue,purple},'markersize',ms,'linewidth',lw);
        %         subplot(3,2,5);
        %         plt.bar(D1.SN,[D1.bFirst,D1.bSecond,D1.b],'style',sty,'leg',sessleg); legend off; title('Sequence recognition test');
        %         xlabel('Subject number'); ylabel('Bias'); xt=xticks; xticks(xt(2:3:end)); xticklabels(unique(D.SN)); set(gca,'fontsize',fs);
        %         subplot(3,2,6);
        %         plt.bar(ones(numel(D1.SN),1),[D1.bFirst,D1.bSecond,D1.b],'style',sty,'leg',sessleg);
        %         xlabel('Group average'); ylabel('Bias'); xticklabels(''); axis square; set(gca,'fontsize',fs); %ylim([0 0.2]);
        
        % out
        varargout={D}; %return main structure
        
    case 'date' % analysis of time difference between day 4 (test) and follow-up session (retest)
        % load excel sheet with demographic information
        xls_fn=fullfile(pathToDropbox,'sr2_logfile.xlsx'); xls_sheet='Demographics';
        [num,text,~]=xlsread(xls_fn,xls_sheet);
        
        % store data in dataframe structure
        idx = logical(num(1:20,5));
        
        D.dateTest = text(2:21, 8);
        D.dateRetest = text(28:47, 8);
        
        D.dateDiffDays = datenum(D.dateRetest(idx), 'dd.mm.yy') - datenum(D.dateTest(idx), 'dd.mm.yy');
        D.dateDiffWeeks = round((datenum(D.dateRetest(idx), 'dd.mm.yy') - datenum(D.dateTest(idx), 'dd.mm.yy')) ./ 7);
        
        % out
        varargout={D}; %return main structure
        
    case 'force' % data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_free_all_data.mat'));
        end
        
        % load force data
        b=0; %initialize block number
        for t=1:numel(D.TN)
            if D.BN(t,1)~=b
                b=D.BN(t,1); %update block number within the loop
                MOV=movload(fullfile(pathToData,sprintf('sr2_s%02d_free_%02d.mov',D.SN(t,1),D.BN(t,1)))); %load MOV file
                fprintf(1,'\nsubj: s%02d   block: %02d\n\n',D.SN(t,1),D.BN(t,1));
            end
            MOV_trial=MOV{D.TN(t,1)};
            time=(2:2:(size(MOV_trial(:,3),1))*2)';
            force=smooth_kernel(MOV_trial(:,4:end),4);
            
            % calculate peak press force
            D.peak_press_force_T(t,1)=max(force(:,RH(1)));
            D.peak_press_force_I(t,1)=max(force(:,RH(2)));
            D.peak_press_force_M(t,1)=max(force(:,RH(3)));
            D.peak_press_force_R(t,1)=max(force(:,RH(4)));
            D.peak_press_force_L(t,1)=max(force(:,RH(5)));
            D.peak_press_force_firstPress(t,1)=max(force(:,RH(D.press1(t,1))));
            D.peak_press_force_avgAllFingers(t,1)=mean(max(force(:,RH)));
            
            % calculate delta time between press and release
            delta_t_pressRelease_T=max(time(force(:,RH(1))>1))-min(time(force(:,RH(1))>1)); if isempty(delta_t_pressRelease_T); D.delta_t_pressRelease_T(t,1)=NaN; else; D.delta_t_pressRelease_T(t,1)=delta_t_pressRelease_T; end
            delta_t_pressRelease_I=max(time(force(:,RH(2))>1))-min(time(force(:,RH(2))>1)); if isempty(delta_t_pressRelease_I); D.delta_t_pressRelease_I(t,1)=NaN; else; D.delta_t_pressRelease_I(t,1)=delta_t_pressRelease_I; end
            delta_t_pressRelease_M=max(time(force(:,RH(3))>1))-min(time(force(:,RH(3))>1)); if isempty(delta_t_pressRelease_M); D.delta_t_pressRelease_M(t,1)=NaN; else; D.delta_t_pressRelease_M(t,1)=delta_t_pressRelease_M; end
            delta_t_pressRelease_R=max(time(force(:,RH(4))>1))-min(time(force(:,RH(4))>1)); if isempty(delta_t_pressRelease_R); D.delta_t_pressRelease_R(t,1)=NaN; else; D.delta_t_pressRelease_R(t,1)=delta_t_pressRelease_R; end
            delta_t_pressRelease_L=max(time(force(:,RH(5))>1))-min(time(force(:,RH(5))>1)); if isempty(delta_t_pressRelease_L); D.delta_t_pressRelease_L(t,1)=NaN; else; D.delta_t_pressRelease_L(t,1)=delta_t_pressRelease_L; end
            delta_t_pressRelease_firstPress=max(time(force(:,RH(D.press1(t,1)))>1))-min(time(force(:,RH(D.press1(t,1)))>1)); if isempty(delta_t_pressRelease_firstPress); D.delta_t_pressRelease_firstPress(t,1)=NaN; else; D.delta_t_pressRelease_firstPress(t,1)=delta_t_pressRelease_firstPress; end
            D.delta_t_pressRelease_avgAllFingers(t,1)=nanmean([D.delta_t_pressRelease_T(t,1),D.delta_t_pressRelease_I(t,1),D.delta_t_pressRelease_M(t,1),D.delta_t_pressRelease_R(t,1),D.delta_t_pressRelease_L(t,1)]);
            
            % calculate baseline force (avg force 3000-0 ms pre-press, first finger only)
            D.baseline_force_T(t,1)=mean(force(1:find(time==3100),RH(1)));
            D.baseline_force_I(t,1)=mean(force(1:find(time==3100),RH(2)));
            D.baseline_force_M(t,1)=mean(force(1:find(time==3100),RH(3)));
            D.baseline_force_R(t,1)=mean(force(1:find(time==3100),RH(4)));
            D.baseline_force_L(t,1)=mean(force(1:find(time==3100),RH(5)));
            D.baseline_force_firstPress(t,1)=mean(force(1:find(time==3100),RH(D.press1(t,1))));
            D.baseline_force_avgAllFingers(t,1)=mean(mean(force(1:find(time==3100),RH)));
            
            % calculate peak press velocity
            if isempty(delta_t_pressRelease_T)||(delta_t_pressRelease_T==0); D.peak_press_velocity_T(t,1)=NaN; else; D.peak_press_velocity_T(t,1)=max(diff(force(force(:,RH(1))>1,RH(1)))); end
            if isempty(delta_t_pressRelease_I)||(delta_t_pressRelease_I==0); D.peak_press_velocity_I(t,1)=NaN; else; D.peak_press_velocity_I(t,1)=max(diff(force(force(:,RH(2))>1,RH(2)))); end
            if isempty(delta_t_pressRelease_M)||(delta_t_pressRelease_M==0); D.peak_press_velocity_M(t,1)=NaN; else; D.peak_press_velocity_M(t,1)=max(diff(force(force(:,RH(3))>1,RH(3)))); end
            if isempty(delta_t_pressRelease_R)||(delta_t_pressRelease_R==0); D.peak_press_velocity_R(t,1)=NaN; else; D.peak_press_velocity_R(t,1)=max(diff(force(force(:,RH(4))>1,RH(4)))); end
            if isempty(delta_t_pressRelease_L)||(delta_t_pressRelease_L==0); D.peak_press_velocity_L(t,1)=NaN; else; D.peak_press_velocity_L(t,1)=max(diff(force(force(:,RH(5))>1,RH(5)))); end
            if isempty(delta_t_pressRelease_firstPress)||(delta_t_pressRelease_firstPress==0); D.peak_press_velocity_firstPress(t,1)=NaN; else; D.peak_press_velocity_firstPress(t,1)=max(diff(force(force(:,RH(D.press1(t,1)))>1,RH(D.press1(t,1))))); end
            D.peak_press_velocity_avgAllFingers(t,1)=nanmean([D.peak_press_velocity_T(t,1),D.peak_press_velocity_I(t,1),D.peak_press_velocity_M(t,1),D.peak_press_velocity_R(t,1),D.peak_press_velocity_L(t,1)]);
            
            % calculate peak release velocity
            if isempty(delta_t_pressRelease_T)||(delta_t_pressRelease_T==0); D.peak_release_velocity_T(t,1)=NaN; else; D.peak_release_velocity_T(t,1)=min(diff(force(force(:,RH(1))>1,RH(1)))); end
            if isempty(delta_t_pressRelease_I)||(delta_t_pressRelease_I==0); D.peak_release_velocity_I(t,1)=NaN; else; D.peak_release_velocity_I(t,1)=min(diff(force(force(:,RH(2))>1,RH(2)))); end
            if isempty(delta_t_pressRelease_M)||(delta_t_pressRelease_M==0); D.peak_release_velocity_M(t,1)=NaN; else; D.peak_release_velocity_M(t,1)=min(diff(force(force(:,RH(3))>1,RH(3)))); end
            if isempty(delta_t_pressRelease_R)||(delta_t_pressRelease_R==0); D.peak_release_velocity_R(t,1)=NaN; else; D.peak_release_velocity_R(t,1)=min(diff(force(force(:,RH(4))>1,RH(4)))); end
            if isempty(delta_t_pressRelease_L)||(delta_t_pressRelease_L==0); D.peak_release_velocity_L(t,1)=NaN; else; D.peak_release_velocity_L(t,1)=min(diff(force(force(:,RH(5))>1,RH(5)))); end
            if isempty(delta_t_pressRelease_firstPress)||(delta_t_pressRelease_firstPress==0); D.peak_release_velocity_firstPress(t,1)=NaN; else; D.peak_release_velocity_firstPress(t,1)=min(diff(force(force(:,RH(D.press1(t,1)))>1,RH(D.press1(t,1))))); end
            D.peak_release_velocity_avgAllFingers(t,1)=nanmean([D.peak_release_velocity_T(t,1),D.peak_release_velocity_I(t,1),D.peak_release_velocity_M(t,1),D.peak_release_velocity_R(t,1),D.peak_release_velocity_L(t,1)]);
            
            % plot trial
            fig=0; % 1=yes|0=no
            if fig==1
                figure;
                if D.hand(t,1)==1
                    plot(time,force(:,LH),'LineWidth',lw);
                    title('Force traces for LEFT hand presses','FontSize',fs);
                elseif D.hand(t,1)==2
                    plot(time,force(:,RH),'LineWidth',lw);
                    title('Force traces for RIGHT hand presses','FontSize',fs);
                end
                drawline(D.pressTime1(t,1),'dir','vert','color',[.8 .8 .8],'linewidth',lw);
                drawline(D.pressTime2(t,1),'dir','vert','color',[.6 .6 .6],'linewidth',lw);
                drawline(D.pressTime3(t,1),'dir','vert','color',[.4 .4 .4],'linewidth',lw);
                drawline(D.pressTime4(t,1),'dir','vert','color',[.2 .2 .2],'linewidth',lw);
                drawline(D.pressTime5(t,1),'dir','vert','color',[.0 .0 .0],'linewidth',lw);
                drawline([3200-D.respWindow(t,1),...
                    3200+D.respWindow(t,1)],'dir','vert','color','r','linewidth',lw);
                xlabel('Time (ms)'); ylabel('Force (N)'); set(gca,'FontSize',fs);
                legend({'Thumb','Index','Middle','Ring','Little','press 1','press 2','press 3','press 4','press 5','resp window'},'FontSize',fs)
            end
        end
        
        % out
        if nargin>1 % single subj data
            save(fullfile(pathToAnalyze,sprintf('sr2_free_force_%s.mat',subj{1})),'-struct', 'D'); %save .mat file
        else % group data
            save(fullfile(pathToAnalyze,sprintf('sr2_free_force_group_data_n=%02d.mat',numel(subj))),'-struct', 'D'); %save .mat file
        end
        varargout={D}; %return main structure
        
    case 'peak_force' % peak press force
        if nargin>1 % load single subj data
            subj=varargin(1);
            s=str2double(subj{1}(2:end));
            D=load(fullfile(pathToAnalyze,sprintf('sr2_free_force_group_data_n=%02d.mat',ns)));
            D=getrow(D,D.SN==s);
        else % load group data
            D=load(fullfile(pathToAnalyze,sprintf('sr2_free_force_group_data_n=%02d.mat',ns)));
        end
        
        %---------------------------------------------------------------------------------------------
        % create summary tables for peak press force
        D2=tapply(D,{'SN','BN','train','free','mix'},...
            {D.peak_press_force_avgAllFingers,'nanmedian','name','F'},...
            'subset',D.isError==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        D2 = normData(D2, {'F'}, 'sub');
        
        %         % open figure
        %         if nargin>1; figure('Name',sprintf('Peak press force - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Peak press force - group (N=%d)',ns)); end
        %         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,2); title('');
        plt.box([D2.mix],D2.normF, 'split',[D2.free D2.train], 'style',tfsty_cbs, 'leg','auto', 'leglocation','northwest');
        xticks([1.75,3.75]); xticklabels({'Blocked', 'Mixed'}); xlabel('Block type'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
        ylim([1.8 5.2]);
        
        %stats
        D2=tapply(D,{'SN','train','free','mix'},...
            {D.peak_press_force_avgAllFingers,'nanmedian','name','F'},...
            'subset',D.isError==0);
        anovaMixed(D2.F, D2.SN, 'within', [D2.mix, D2.free, D2.train], {'mix', 'free', 'train'});
        ttest(D2.F(D2.mix==0 & D2.free==0 & D2.train==0), D2.F(D2.mix==0 & D2.free==0 & D2.train==1), 2, 'paired');
        ttest(D2.F(D2.mix==0 & D2.free==1 & D2.train==0), D2.F(D2.mix==0 & D2.free==1 & D2.train==1), 2, 'paired');
        ttest(D2.F(D2.mix==1 & D2.free==0 & D2.train==0), D2.F(D2.mix==1 & D2.free==0 & D2.train==1), 2, 'paired');
        ttest(D2.F(D2.mix==1 & D2.free==1 & D2.train==0), D2.F(D2.mix==1 & D2.free==1 & D2.train==1), 2, 'paired');
%         %---------------------------------------------------------------------------------------------
%         % create summary tables for peak press force
%         D2=tapply(D,{'SN','BN','train'},...
%             {D.peak_press_force_avgAllFingers,'nanmedian','name','F'},...
%             'subset',D.isError==0 & D.free==0);
%         
%         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
%         D2 = normData(D2, {'F'}, 'sub');
%         
%         %         % open figure
%         %         if nargin>1; figure('Name',sprintf('Peak press force - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Peak press force - group (N=%d)',ns)); end
%         %         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
%         
%         subplot(2,2,3); title('Forced-RT');
%         plt.box(D2.BN,D2.normF, 'split',D2.train, 'style',trsty_cbs, 'leg',trleg, 'leglocation','northwest');
%         xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
%         ylim([1.8 5.2]);
%         
%         %---------------------------------------------------------------------------------------------
%         % create summary tables for peak press force
%         D2=tapply(D,{'SN','BN','train'},...
%             {D.peak_press_force_avgAllFingers,'nanmedian','name','F'},...
%             'subset',D.isError==0 & D.free==1);
%         
%         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
%         D2 = normData(D2, {'F'}, 'sub');
%         
%         subplot(2,2,4); title('Free-RT');
%         plt.box(D2.BN,D2.normF, 'split',D2.train, 'style',trsty_cbs, 'leg',trleg, 'leglocation','northwest');
%         xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
%         ylim([1.8 5.2]);
        
        % out
        varargout={D}; %return main structure
        
    case 'press_dur_ms' % analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
        if nargin > 1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr2_%s.mat',subj)) ); %load data for this subject
            D.SN = ones( numel(D.TN), 1) * str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr2_free_all_data.mat') );
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Press duration / overlap - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Press duration / overlap - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % create summary table
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {sum(D.press2rel_dur,2),'nanmedian', 'name','p2r'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T1.p2r           = reshape(T.p2r, [], 1);
        T1.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T1.BN            = repmat(T.BN, numel(T.p2r(1,:)), 1);
        T1.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T1.day           = repmat(T.day, numel(T.p2r(1,:)), 1);
        T1.P             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T1 = normData(T1, {'p2r'}, 'sub');
        
        subplot(2,2,1); title(''); hold on;
        plt.line([T1.day T1.BN], T1.normp2r, 'plotfcn','nanmean', 'split',T1.train, 'errorbars','shade', 'style',dursty, 'leg',trleg, 'leglocation','northeast');
        xlabel('Block number'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        
        % create summary table
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {sum(D.rel2press_dur,2),'nanmedian', 'name','r2p'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T3.r2p           = reshape(T.r2p, [], 1);
        T3.SN            = repmat(T.SN, numel(T.r2p(1,:)), 1);
        T3.BN            = repmat(T.BN, numel(T.r2p(1,:)), 1);
        T3.train         = repmat(T.train, numel(T.r2p(1,:)), 1);
        T3.day           = repmat(T.day, numel(T.r2p(1,:)), 1);
        T3.T             = reshape(repmat(1:numel(T.r2p(1,:)), numel(T.r2p(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T3 = normData(T3, {'r2p'}, 'sub');
        
        subplot(2,2,1); title(''); hold on;
        [x,~,~]=plt.line([T3.day T3.BN], T3.normr2p, 'plotfcn','nanmean', 'split',T3.train, 'errorbars','shade', 'style',delsty, 'leg',trleg, 'leglocation','northeast');
        xlabel('Block number'); ylabel('Duration (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        
        hold on;
        subplot(2,2,1);
        drawline(x([1 12, 13 24, 25 36, 37 48]),'dir','vert','linestyle',':','color','k');
        drawline(0, 'dir','horz', 'linestyle',':', 'color','k');
        hold off;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end
end