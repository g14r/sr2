function[varargout]=sr2_training_analyze(what,varargin)
%% [varargout]=sr2_training_analyze(what,varargin)
% SequenceRepetition experiment (sr2_training) analysis of behavioral data
%
% usage|example calls:
%
%                       sr2_training_analyze('all_subj');                               %pre-analysis: run the subject routine for all_subj
%                       sr2_training_analyze('all_subj',{'s01'});                       %pre-analysis: run the subject routine for selected subjects
%                       [all_data]=sr2_training_analyze('all_data');                    %pre-analysis: create .mat file containing data from subjects in subj
%
%                       [D]=sr2_training_analyze('rtt');                                %group          analysis of reaction time task (RT), for different prep times
%                       [D,rtt_x_day]=sr2_training_analyze('rtt','s01');                %single subject analysis of reaction time task (RT), for different prep times
%                       [D]=sr2_training_analyze('rtt_logfit');                         %single subject results of logistic fit estimation (given a certain ACC/y, the corresponding RT/x)
%
%                       [D]=sr2_training_analyze('rtt_rpt');                            %group          analysis of reaction time task (rtt), with real preparation time (rpt)
%                       [D,rtt_x_day]=sr2_training_analyze('rtt_rpt','s01');            %single subject analysis of reaction time task (rtt), with real preparation time (rpt)
%                       [D]=sr2_training_analyze('rtt_rpt_logfit');                     %single subject analysis of reaction time task (rtt), with real preparation time (rpt), per subject (subj), including logfit
%
%                       [D]=sr2_training_analyze('training');                           %group          analysis of training effects over time on MT and ER, for different prep times
%                       [D]=sr2_training_analyze('training','s01');                     %single subject analysis of training effects over time on MT and ER, for different prep times
%
%                       [D]=sr2_training_analyze('prepTime_works');                     %group          analysis of the effects of prep time on MT and IPI, early in training
%                       [D]=sr2_training_analyze('prepTime_works','s01');               %single subject analysis of the effects of prep time on MT and IPI, early in training
%
%                       [D]=sr2_training_analyze('prepTime');                           %group          analysis of the effects of prep time on MT and ER, for traingved/random sequences
%                       [D]=sr2_training_analyze('prepTime','s01');                     %single subject analysis of the effects of prep time on MT and ER, for trained/random sequences
%
%                       [D]=sr2_training_analyze('points');                             %group          analysis of the proportion of points/errors over time (effect of thresholding)
%                       [D]=sr2_training_analyze('points','s01');                       %single subject analysis of the proportion of points/errors over time (effect of thresholding)
%
%                       [D]=sr2_training_analyze('RT');                                 %group          analysis of the distribution of reaction times
%                       [D]=sr2_training_analyze('RT','s01');                           %single subject analysis of the distribution of reaction times
%
%                       [D]=sr2_training_analyze('seq');                                %group          analysis of the effects of different sequence types on behavioral measures
%                       [D]=sr2_training_analyze('seq','s01');                          %single subject analysis of the effects of different sequence types on behavioral measures
%
%                       [D]=sr2_training_analyze('IPI');                                %group          analysis of the inter-press intervals (IPIs) as a function of training and prep time
%                       [D]=sr2_training_analyze('IPI','s01');                          %single subject analysis of the inter-press intervals (IPIs) as a function of training and prep time
%
%                       [D]=sr2_training_analyze('force');                              %group          data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
%                       [D]=sr2_training_analyze('force','s01');                        %single subject data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
%
%                       [D]=sr2_training_analyze('peak_force');                         %group          peak press force
%                       [D]=sr2_training_analyze('peak_force','s01');                   %single subject peak press force
%
%                       [D]=sr2_training_analyze('press_time');                         %group          delta time press-release
%                       [D]=sr2_training_analyze('press_time','s01');                   %single subject delta time press-release
%
%                       [D]=sr2_training_analyze('baseline_force');                     %group          pre-press baseline force
%                       [D]=sr2_training_analyze('baseline_force','s01');               %single subject pre-press baseline force
%
%                       [D]=sr2_training_analyze('peak_velocity');                      %group          peak press and release velocity
%                       [D]=sr2_training_analyze('peak_velocity','s01');                %single subject peak press and release velocity
%
%                       [D]=sr2_training_analyze('press_dur_ET');                       %group          analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur_ET','s01');                 %single subject analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur_ms');                       %group          analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur_ms','s01');                 %single subject analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur_PT');                       %group          analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur_PT','s01');                 %single subject analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%
% --
% gariani@uwo.ca - 2018.01.22

%% globals

% paths
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr2';
pathToAnalyze='/Users/gariani/Documents/data/SequenceRepetition/sr2/analyze';

% subjects
subj={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20'};
ns=numel(subj);
subvec=zeros(1,ns);
for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
darkblue=[0,90,140]/255;
blue=[49,130,189]/255;
lightblue=[145,195,225]/255;
red=[222,45,38]/255;
lightred=[252,146,114]/255;
green=[49,163,84]/255;
lightgreen=[161,217,155]/255;
darkorange=[210,80,0]/255;
orange=[253,141,60]/255;
lightorange=[250,170,100]/255;
darkyellow=[230,190,70]/255;
yellow=[254,196,79]/255;
lightyellow=[255,237,160]/255;
purple=[117,107,177]/255;
lightpurple=[188,189,220]/255;
darkgray = [80,80,80]/255;
gray3 = [50,50,50]/255;
gray2 = [100,100,100]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
ptc_un_cbs={[229	158	102]/255;[221	126	50]/255;[213	94	0]/255;[149	65	0]/255};
ptc_tr_cbs={[153	198	224]/255;[76	156	201]/255;[0	102	160]/255;[0	68	106]/255};
fingc={green,yellow,red,purple,blue};
trc={lightgray,black};

% plot defaults
fs=20; %default font size for all figures
lw=4; %2; %2;%3; %4; %3; %default line width for all figures
ms=12; %6; %8;%12; %10; %default marker size for all figures

% styles
style.reset;
style.custom({darkblue,blue,lightblue,red,lightred,darkorange,orange,lightorange,darkyellow,yellow,lightyellow,purple,lightpurple,gray,gray2,gray3,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
trsty_cbs=style.custom({cbs_red,cbs_blue},'markersize',ms,'linewidth',lw,'errorwidth',lw);
unday1sty_cbs=style.custom({cbs_red,cbs_blue},'markersize',ms,'linewidth',lw,'linestyle',':','markersize',ms,'errorwidth',lw);%,'markertype','none');
ptUNsty=style.custom(ptc_un_cbs,'markersize',ms,'linewidth',lw);
ptTRsty=style.custom(ptc_tr_cbs,'markersize',ms,'linewidth',lw);
ipisty=style.custom({lightgray,gray,darkgray,black},'markersize',ms,'linewidth',lw,'errorwidth',lw);
delsty=style.custom({lightorange,lightblue},'markersize',ms,'linewidth',lw);
dursty=style.custom({darkorange,darkblue},'markersize',ms,'linewidth',lw);
rttModelsty=style.custom({lightgray,gray,darkgray,black},'linewidth',lw, 'markertype','none', 'errorbars','shade');

% legends
trleg={'Untrained','Trained'};
ptleg={'400','800','1600','2400'};
fingleg={'Thumb','Index','Middle','Ring','Little'};
dleg={'Day 1','Day 2','Day 3','Day 4'};

% hands
LH=[1 2 3 4 5]; %left hand column indices (force analysis)
RH=[6 7 8 9 10]; %right hand column indices (force analysis)

%% types of analysis
switch (what)
    case 'all_subj' % pre-analysis: run the subject routine for all_subj
        if nargin>1; subj=varargin{1}; end
        for s=1:numel(subj)
            sr2_training_subj(subj{s},0); % run sr2_training_subj.m routine (without plot)
        end
        
    case 'all_data' % pre-analysis: create .mat file containing data from all subjects
        all_data=[];
        for s=1:ns
            fprintf('\n%s\n\n',subj{s});
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{s}))); %load data structure for each subject
            D.SN=ones(numel(D.TN),1)*subvec(s); % add information about subject number
            % add BN info per session
            for iDay=1:numel(unique(D.day))
                nBlocksPerDay=numel(unique(D.BN(D.day==iDay))); %all blocks including trained/untrained
                D.BN_day(D.day==iDay,1)=D.BN(D.day==iDay,1)-(nBlocksPerDay*(iDay-1));
            end
            %D.BN_session2=(ceil((1:numel(D.BN))/(numel(D.BN)/nBlocksOverall)))';
            all_data=addstruct(all_data,D); % append data structures from each subject
        end
        save(fullfile(pathToAnalyze,'sr2_training_all_data.mat'),'-struct', 'all_data'); %save all_data.mat file
        varargout={all_data};
        
    case 'rtt' % analysis of reaction time task (RT), for different prep times
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % ----------------------------------------------------------------------------------------------------------------------------
        % create summary tables for ACC
        T = tapply(D,{'SN','prepTime','day'},...
            {(1-D.pressError)*100,'nanmean','name','ACC'}, ...
            'subset',D.dummy==0 & D.rtt==1 & D.timingError==0);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.prepTime, T.day, T.normACC, 'length');
        
        % open figure
        if nargin>1; figure('Name',sprintf('RTT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('RTT - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1);
        %[x,y]=plt.line(T.prepTime,T.normACC, 'errorbars','shade','split',T.day,'style',rttDatasty,'leg',dleg,'leglocation','east');
        [x,y]=plt.line(T.prepTime,T.normACC, 'errorbars','shade','split',T.day,'style',ipisty,'leg',dleg,'leglocation','east');
        title('Data');
        xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
        xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
        axis square; set(gca,'fontsize',fs);
        %chance
        drawline(20,'dir','horz','linestyle','-','color',black);
        %ref points
        drawline(400,'dir','vert','linestyle',':','color',black);
        drawline(y(1,5),'dir','horz','linestyle',':','color',black);
        drawline(y(4,5),'dir','horz','linestyle',':','color',black);
        hold off;
        
        %----------------------------------------------------------------------------------------------------------------------------
        % log fit
        rtt_x_day=zeros(numel(unique(T.day)),1);
        rtt_y=0.80; % the accuracy value (from 0 to 1) that you want to use to determine subj RT with that accuracy (from logistic fit)
        % fit logistic function
        R.x=[]; R.y_hat=[]; R.day=[]; y=y/100;
        for iDay=1:numel(unique(T.day))
            %theta_zero=[0.1,400,.20]';
            theta_zero=[0.01, 5]';
            [theta_hat]=fitlog(x,y(iDay,:),theta_zero);
            a=theta_hat(1); b=theta_hat(2);
            %c=theta_hat(3);
            c=.20;
            rtt_x_day(iDay)= (b - log((1-c)/(rtt_y-c) - 1)) / a;
            y_hat=modlog(theta_hat,x);
            R.x=[R.x;x];
            R.y_hat=[R.y_hat;y_hat];
            R.day=[R.day;ones(numel(unique(T.prepTime)),1)*iDay];
        end
        subplot(2,2,2);
        plt.line(R.x,R.y_hat*100,'split',R.day,'style',ipisty,'leg',dleg,'leglocation','east');
        title('Logistic fit');
        xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
        xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
        axis square; set(gca,'fontsize',fs);
        %chance
        drawline(20,'dir','horz','linestyle','-','color',black);
        %ref points
        drawline(rtt_y*100,'dir','horz','linestyle',':','color',black);
        drawline(rtt_x_day(1),'dir','vert','linestyle',':','color',black);
        drawline(rtt_x_day(4),'dir','vert','linestyle',':','color',black);
        hold off;
        
        % stats
        ST = tapply(D,{'SN', 'day'},...
            {(1-D.pressError)*100,'nanmean','name','ACC', 'subset',D.prepTime==400}, ...
            'subset',D.dummy==0 & D.rtt==1 & D.timingError==0);
        ttest(ST.ACC(ST.day==1), ST.ACC(ST.day==4), 2, 'paired');
        
        % out
        D.T=T; D.rtt_y=rtt_y; D.rtt_x_day=rtt_x_day; D.R=R; %incorporate the sub-structures as fields of main structure
        varargout={D,rtt_x_day}; %return main structure
        
    case 'rtt_logfit' % single subject results of logistic fit estimation (given a certain ACC/y, the corresponding RT/x)
        subj={'s01','s02','s03','s04','s05','s06',      's08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20'}; % subj 7 does not have enough trials in short prep times for model fitting
        ns=numel(subj);
        subvec=zeros(1,ns); for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end
        if ~exist(fullfile(pathToAnalyze, 'sr2_model_fits_subj.mat'), 'file')
            D.model_fits_x = []; D.model_fits_y = []; D.model_fits_day = []; D.model_fits_SN = []; rtt_x_day = zeros(ns, 4);
            for s=1:ns
                sn=subj{s};
                S=sr2_training_analyze('rtt',sn); close all;
                rtt_x_day(s,:)=real(S.rtt_x_day);
                % store info about model fits per subject
                D.model_fits_x      = [D.model_fits_x; S.R.x];
                D.model_fits_y      = [D.model_fits_y; S.R.y_hat];
                D.model_fits_day    = [D.model_fits_day; S.R.day];
                D.model_fits_SN     = [D.model_fits_SN; ones(numel(S.R.x),1) * subvec(s)];
            end
            D.rtt_x_all = rtt_x_day;
            D.rtt_x = mean(rtt_x_day, 1);
            D.rtt_y = S.rtt_y;
            save( fullfile(pathToAnalyze, 'sr2_model_fits_subj.mat'), '-struct', 'D'); %save all_data.mat file
        else
            [D] = load( fullfile(pathToAnalyze, 'sr2_model_fits_subj.mat') );
        end
        
        %---------------------------------------------------------------------------------------------------
        % produce old plot
        sr2_training_analyze('rtt'); hold on;
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        D.SN = D.model_fits_SN;
        D = normData(D, {'model_fits_y'}, 'sub');
        
        subplot(2,2,3);
        %[~,~] = plt.line(D.model_fits_x, D.model_fits_y*100, 'errorbars','shade', 'split',D.model_fits_day, 'style',rttModelsty, 'leg',dleg, 'leglocation','east');
        [~,~] = plt.line(D.model_fits_x, D.normmodel_fits_y*100, 'errorbars','shade', 'split',D.model_fits_day, 'style',ipisty, 'leg',dleg, 'leglocation','east');
        title('Logistic fit'); hold on;
        xticklabels({'200','','300','','400','','500','','600',''}); xlim([180 670]); ylim([11 100]);
        xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
        axis square; set(gca,'fontsize',fs);
        %chance
        drawline(20,'dir','horz','linestyle','-','color',black);
        
        %ref points
        % gain in planning speed
        drawline(D.rtt_y*100, 'dir','horz', 'linestyle',':', 'color',black);
        drawline(D.rtt_x(1), 'dir','vert', 'linestyle',':', 'color',black);
        drawline(D.rtt_x(4), 'dir','vert', 'linestyle',':', 'color',black);
        hold off;
        
        %---------------------------------------------------------------------------------------------------
        % stats
        % gain in selection accuracy
        ST = tapply(D,{'model_fits_SN', 'model_fits_day'},...
            {D.model_fits_y*100, 'nanmean', 'name','y', 'subset',D.model_fits_x==400});
        ttest(ST.y(ST.model_fits_day==1), ST.y(ST.model_fits_day==4), 2, 'paired');
        
        % gain in planning speed
        ttest(D.rtt_x_all(:, 1), D.rtt_x_all(:, 4), 2, 'paired');
        
        % out
        varargout={D};
        
    case 'training' % analysis of training effects over time on MT and ER, for different prep times
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        %---------------------------------------------------------------------------------------------------
        % create summary tables for MT
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {D.MT,'nanmedian', 'name','MT'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Training - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Training - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.MT, 'length');
        
        subplot(2,2,1); title('Learning - SET');
        [x,~] = plt.line([T.day T.BN],T.normMT, 'split',T.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Block number (overall)'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        hold on; drawline(x([1 12, 13 24, 25 36, 37 48]),'dir','vert','linestyle',':','color','k'); hold off;
        
        % stats
        % train vs untrain, day 4
        T = tapply(D,{'SN', 'train'},...
            {D.MT,'mean', 'name','MT', 'subset', D.day==4 & D.BN_day>10},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        ttest(T.MT(T.train==0), T.MT(T.train==1), 2, 'paired');
        [T.t, T.p] = ttest(T.MT(T.train==0), T.MT(T.train==1), 2, 'paired');
        
        % day 1 vs day 4, untrain
        T = tapply(D,{'SN','day'},...
            {D.MT,'mean', 'name','MT', 'subset', D.train==0 & (D.day==1 | D.day==4)},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        ttest(T.MT(T.day==1), T.MT(T.day==4), 2, 'paired');
        [T.t, T.p] = ttest(T.MT(T.day==1), T.MT(T.day==4), 2, 'paired');
        
        %---------------------------------------------------------------------------------------------------
        % create summary tables for ACC
        T = tapply(D,{'SN','day','train'},...
            {(1-D.isError)*100,'nanmean','name','ACC'}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.BN_day>10);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'ACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.ACC, 'length');
        
        subplot(2,2,2); hold on; title('Learning - ACC');
        plt.line(T.day,T.normACC, 'errorbars','shade', 'split',T.train, 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Training day'); ylabel('Accuracy (%)'); set(gca,'fontsize',fs); xlim([.5 4.5]); ylim([60 100]); axis square;
        
        %---------------------------------------------------------------------------------------------------
        % create summary tables for pressACC
        T = tapply(D,{'SN','day','train'},...
            {(1-D.pressError)*100,'nanmean','name','pACC', 'subset',D.timingError==0}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.BN_day>10);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'pACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.pACC, 'length');
        
        subplot(2,2,3); hold on; title('Learning - press ACC');
        plt.line(T.day,T.normpACC, 'errorbars','shade', 'split',T.train, 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Training day'); ylabel('Press accuracy (%)'); set(gca,'fontsize',fs); xlim([.5 4.5]); ylim([60 100]); axis square;
        
        %---------------------------------------------------------------------------------------------------
        % create summary tables for timingACC
        T = tapply(D,{'SN','day','train'},...
            {(1-D.timingError)*100,'nanmean','name','tACC'}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.BN_day>10);
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'tACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.tACC, 'length');
        
        subplot(2,2,4); hold on; title('Learning - timing ACC');
        plt.line(T.day,T.normtACC, 'errorbars','shade', 'split',T.train, 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Training day'); ylabel('Timing accuracy (%)'); set(gca,'fontsize',fs); xlim([.5 4.5]); ylim([60 100]); axis square;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'prepTime_works' % analysis of the effects of prep time on MT and IPI, early in training
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % open figure - prep time works
        if nargin > 1
            figure('Name',sprintf('Prep time works - subj %02d', str2double(varargin{1}(2:3))));
        else
            figure('Name',sprintf('Prep time works - group (N=%d)', ns));
        end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %---------------------------------------------------------------------------------------------------
        % create summary table for MT
        T = tapply(D, {'SN','prepTime','day'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.SN, T.prepTime, T.MT, 'mean');
        
        %---------------------------------------------------------------------------------------------------
        % plot data
        subplot(2,2,1); title('');
        %plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'errorbars','plusminus_wocap', 'style',darkgraysty);
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'errorbars','plusminus_wocap', 'split',T.day, 'style',ipisty, 'leg',dleg);
        xlabel('Preparation time (ms)'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([1000 1700]);
        hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        %---------------------------------------------------------------------------------------------------
        % stats
        anovaMixed(T.MT, T.SN, 'within', [T.prepTime, T.day], {'prepTime','day'});
        T2 = tapply(D, {'SN','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.BN_day>10);
        anovaMixed(T2.MT, T2.SN, 'within', [T2.prepTime], {'prepTime'});
        %         ttest(T.MT(T.prepTime==400), T.MT(T.prepTime==800), 2, 'paired');
        %         ttest(T.MT(T.prepTime==800), T.MT(T.prepTime==1600), 2, 'paired');
        %         ttest(T.MT(T.prepTime==1600), T.MT(T.prepTime==2400), 2, 'paired');
        
        %---------------------------------------------------------------------------------------------------
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        
        % create summary table for MT
        T=tapply(D,{'SN','prepTime'},...
            {D.IPI_1,'nanmean','name','IPI1'},...
            {D.IPI_2,'nanmean','name','IPI2'},...
            {D.IPI_3,'nanmean','name','IPI3'},...
            {D.IPI_4,'nanmean','name','IPI4'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.BN_day>10);
        for i=1:size(D.IPI,2)
            T.IPI(:,i)=eval(sprintf('T.IPI%d',i));
            T=rmfield(T,sprintf('IPI%d',i));
            T.IPInum(:,i)=ones(size(T.SN,1),1)*i;
            T.SN(:,i)=T.SN(:,1);
            T.prepTime(:,i)=T.prepTime(:,1);
        end
        T.IPI=reshape(T.IPI,size(T.IPI,1)*size(T.IPI,2),1);
        T.IPInum=reshape(T.IPInum,size(T.IPInum,1)*size(T.IPInum,2),1);
        T.SN=reshape(T.SN,size(T.SN,1)*size(T.SN,2),1);
        T.prepTime=reshape(T.prepTime,size(T.prepTime,1)*size(T.prepTime,2),1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.IPInum, T.prepTime, T.normIPI, 'length');
        
        subplot(2,2,2); title('Average all days');
        plt.line(T.IPInum,T.normIPI, 'split',T.prepTime, 'errorbars','shade', 'style',ipisty, 'leg',ptleg, 'leglocation','northeast');
        xlabel('IPI number'); ylabel('IPI duration (ms)'); set(gca,'fontsize',fs); axis square; ylim([130 470]);
        
        %---------------------------------------------------------------------------------------------------
        % stats
        ttest(T.IPI(T.IPInum==1 & T.prepTime==400), T.IPI(T.IPInum==1 & T.prepTime==800), 2, 'paired');
        ttest(T.IPI(T.IPInum==2 & T.prepTime==400), T.IPI(T.IPInum==2 & T.prepTime==800), 2, 'paired');
        ttest(T.IPI(T.IPInum==3 & T.prepTime==400), T.IPI(T.IPInum==3 & T.prepTime==2400), 2, 'paired');
        ttest(T.IPI(T.IPInum==4 & T.prepTime==400), T.IPI(T.IPInum==4 & T.prepTime==800), 2, 'paired');
        
        % out
        D.T = T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'prepTime' % analysis of the effects of prep time on MT and ER, for trained/random sequences
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % open figure
        if nargin > 1
            figure('Name',sprintf('Prep time - subj %02d', str2double(varargin{1}(2:3))));
        else
            figure('Name',sprintf('Prep time - group (N=%d)', ns));
        end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %---------------------------------------------------------------------------------------------------
        % create summary table for MT (day 4)
        %T = tapply(D, {'SN','train','prepTime','seqNum'}, ...
        T = tapply(D, {'SN','train','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        % compute the percentage of MT at prep time 2400 to show the interaction
        %h=figure; [~,y]=plt.xy(T.realPrepTime,T.normMT,T.prepTime,'split',T.train); close(h);
        
        h=figure; [~,y]=plt.xy(T.realPrepTime,(1./T.normMT),T.prepTime,'split',T.train); close(h);
        %h=figure; [~,y1]=plt.xy(T.realPrepTime,(1./T.normMT),T.prepTime,'subset',T.train==0); close(h);
        %h=figure; [~,y2]=plt.xy(T.realPrepTime,(1./T.normMT),T.prepTime,'split',T.seqNum,'subset',T.train==1); close(h);
        
        T.pctMT = zeros(numel(T.normMT),1);
        tridx = T.train==1;
        %tridx = unique(T.seqNum(T.train==1));
        unidx = T.train==0;
        T.pctMT(unidx) = (1 ./ (T.normMT(unidx)))  ./  y(1,end) * 100;
        T.pctMT(tridx) = (1 ./ (T.normMT(tridx)))  ./  y(2,end) * 100;
        
        %---------------------------------------------------------------------------------------------------
        % day 1
        T2 = tapply(D, {'SN','train','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==1 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'MT'}, 'sub');
        
        % compute the percentage of MT at prep time 2400 to show the interaction
        h=figure; [~,y]=plt.xy(T2.realPrepTime,(1./T2.normMT),T2.prepTime,'split',T2.train); close(h);
        T2.pctMT = zeros(numel(T2.normMT),1);
        tridx = T2.train==1;
        unidx = T2.train==0;
        T2.pctMT(unidx) = (1 ./ (T2.normMT(unidx)))  ./  y(1,4) * 100;
        T2.pctMT(tridx) = (1 ./ (T2.normMT(tridx)))  ./  y(2,4) * 100;
        
        %---------------------------------------------------------------------------------------------------
        % plot data
        subplot(2,2,1); title('Day 4 (end of training)');
        plt.xy(T2.realPrepTime,T2.normMT,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs, 'subset',T2.train==0, 'leg','skip');
        %plt.xy(T2.realPrepTime,T2.normMT,T2.prepTime, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs);
        hold on;
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([850 1850]);
        hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        subplot(2,2,2); title('Interaction');
        plt.xy(T.realPrepTime,T.pctMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','east');
        hold on;
        plt.xy(T2.realPrepTime,T2.pctMT,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs, 'subset',T2.train==0, 'leg','skip');
        %plt.xy(T2.realPrepTime,T2.pctMT,T2.prepTime, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs);
        hold on;
        xlabel('Preparation time (ms)'); ylabel('% of ET at 2400'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(tridx))); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([75 105]);
        hold on; drawline(unique(T.prepTime(tridx)),'dir','vert','linestyle',':','color','k'); drawline(100,'dir','horz','linestyle','--','color','k'); hold off;
        
        %---------------------------------------------------------------------------------------------------
        % stats
        T.res = anovaMixed(T.MT,T.SN,'within', [T.prepTime,T.train], {'prepTime','train'});
        ttest(T.MT(T.prepTime==2400 & T.train==0), T.MT(T.prepTime==2400 & T.train==1), 2, 'paired');
        %[T.t, T.p] = ttest(T.MT(T.prepTime==2400 & T.train==0), T.MT(T.prepTime==2400 & T.train==1), 2, 'paired');
        
        T3 = tapply(D, {'SN','prepTime','day'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.train==0 & ismember(D.day,[1,4]));
        T3.res = anovaMixed(T3.MT,T3.SN,'within', [T3.prepTime,T3.day], {'prepTime','day'});
        
        % out
        D.T = T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'points' % analysis of the proportion of points/errors over time (effect of thresholding)
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        % define colors/legend
        c={red,lightred,blue,lightblue,green,lightgreen};
        l={'-3: large timing error','-1: small timing error','0: press error','+1: below avg speed','+2: avg speed','+3: above avg speed'};
        %         c={red,blue,green};
        %         l={'-1: timing error','0: press error','+1: correct'};
        %         D.points(D.points<0)=-1; D.points(D.points>0)=1;
        D1=tapply(D,{'SN','BN','day','points'},...
            {D.points,'numel','name','nPoints'},...
            'subset',D.dummy==0 & D.rtt==0 & D.train==1);
        % open figure
        if nargin>1; figure('Name',sprintf('Points over blocks - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Points over blocks - group (N=%d)',ns)); end
        lineplot([D1.day D1.BN],D1.nPoints,'split',D1.points,'errorbars','shade','shadecolor',c,'linewidth',lw,'linecolor',c,'leg',l,'leglocation','northeast');
        title('Points analysis'); xlabel('Block number'); ylabel('N trials'); set(gca,'fontsize',fs);
        % open figure
        if nargin>1; figure('Name',sprintf('Proportion of points - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Proportion of points - group (N=%d)',ns)); end
        sub=zeros(length(subj),numel(c)); slabels=cell(length(subj),1);
        for s=1:ns
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{s}))); %D.points(D.points<0)=-1; D.points(D.points>0)=1;
            D=tapply(D,{'points'},...
                {D.points,'numel','name','nPoints'},...
                'subset',D.dummy==0 & D.rtt==0 & D.train==1);
            subplot(2,ceil(ns/2),s);
            barplot((1:numel(unique(D.points)))',D.nPoints/sum(D.nPoints),'split',D.points,'facecolor',c); xticklabels({'-3','-1','0','+1','+2','+3'});
            title(sprintf('Subject %02d',subvec(s))); xlabel('Points'); ylabel('Proportion'); set(gca,'fontsize',fs); set(gca,'ylim',[0 1]);
            sub(s,:)=(D.nPoints/sum(D.nPoints))'; slabels{s,:}=sprintf('subj %02d',subvec(s));
        end
        % open figure
        if nargin>1; figure('Name',sprintf('Proportion of points - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Proportion of points - group (N=%d)',ns)); end
        bar(sub,'stacked'); colorbar('Ticks',1:6,'TickLabels',l); map=zeros(numel(c),3); for ic=1:numel(c); map(ic,:)=c{ic}; end; colormap(map);
        title('Proportion of points'); xlabel(''); ylabel('Proportion of trials'); xlim([0.5 ns+0.5]); xticks(1:ns); xticklabels(slabels); set(gca,'fontsize',fs);
        % out
        D.D1=D1; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'RT' % analysis of the distribution of reaction times
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        % RTT
        x=500; %define x-axis range/cut-off
        D1=getrow(D,D.dummy==0 & D.rtt==1 & D.RT<=x & D.RT>=-x); ptimes=[200,300,400,500,600];
        for ipt=1:numel(ptimes); D1.prepTime(D1.prepTime==ptimes(ipt) | D1.prepTime==(ptimes(ipt)+50))=ptimes(ipt); end
        if nargin>1; figure('Name',sprintf('RTT RT-prepTime histogram - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('RTT RT-prepTime histogram - group (N=%d)',ns)); end
        for ipt=1:numel(ptimes)
            subplot(numel(ptimes),1,ipt);
            nurt=numel(unique(D1.RT(D1.prepTime==ptimes(ipt))));
            title(sprintf('Prep time: %d',ptimes(ipt))); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]);
            plt.hist(D1.RT,'numcat',nurt,'percent',1,'subset',D1.prepTime==ptimes(ipt)); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        end
        % STRST
        x=500; %define x-axis range/cut-off
        D2=getrow(D,D.dummy==0 & D.rtt==0 & D.RT<=x & D.RT>=-x); ptimes=unique(D2.prepTime);
        if nargin>1; figure('Name',sprintf('STRST RT-prepTime histogram - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('STRST RT-prepTime histogram - group (N=%d)',ns)); end
        for ipt=1:numel(ptimes)
            subplot(numel(ptimes),1,ipt);
            nurt=numel(unique(D2.RT(D2.prepTime==ptimes(ipt))));
            title(sprintf('Prep time: %d',ptimes(ipt))); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]);
            plt.hist(D2.RT,'numcat',nurt,'percent',1,'subset',D2.prepTime==ptimes(ipt)); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        end
        % STRST individual subjects, collapsing prepTime
        x=500; %define x-axis range/cut-off
        D3=getrow(D,D.dummy==0 & D.rtt==0 & D.RT<=x & D.RT>=-x);
        if nargin>1; figure('Name',sprintf('STRST RT histogram - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('STRST RT histogram - group (N=%d)',ns)); end
        for is=1:ns
            subplot(ceil(ns/2),2,is);
            nurt=numel(unique(D3.RT(D3.SN==subvec(is))));
            title(sprintf('Subject %02d',subvec(is))); xlabel('RT (ms)'); ylabel('Trials (%)'); set(gca,'fontsize',fs); xlim([-x x]);
            plt.hist(D3.RT,'numcat',nurt,'percent',1,'subset',D3.SN==subvec(is)); drawline([-100,100],'dir','vert','linestyle','--','color','r','linewidth',2);
        end
        %         % Scatter plot MT-RT
        %         if nargin>1; figure('Name',sprintf('Scatterplot MT-RT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Scatterplot RT-MT - group (N=%d)',ns)); end
        %         scatterplot(D.MT,D.RT,'split',D.train,'regression','linear','markertype','o','markerfill','auto','markersize',9,'printcorr','subset',D.dummy==0 & D.isError==0 & D.rtt==0 & D.RT<=1000 & D.MT<=3000);
        %         title('Relationship RT-MT'); xlabel('MT (ms)'); ylabel('RT (ms)'); set(gca,'fontsize',fs); %ylim([-500 1000]);
        % out
        D.D1=D1; D.D2=D2; D.D3=D3; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'seq' % analysis of the distribution of reaction times
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        % create summary tables for MT and press ER
        D1=tapply(D,{'SN','seqNum'},...
            {D.MT,'mean','name','MT','subset',D.isError==0},...
            {D.pressError,'mean','name','pER','subset',D.timingError==0},...
            'subset',D.dummy==0 & D.rtt==0 & D.train==1);
        % create summary tables for MT and press ER
        D2=tapply(D,{'SN','day','seqNum'},...
            {D.MT,'mean','name','MT','subset',D.isError==0},...
            {D.pressError,'mean','name','pER','subset',D.timingError==0},...
            'subset',D.dummy==0 & D.rtt==0 & D.train==1);
        seqlabels={'1A','1B','2A','2B','3A','3B','4A','4B','5A','5B'};
        % MT
        if nargin>1; figure('Name',sprintf('Seq MT - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Seq MT - group (N=%d)',ns)); end
        subplot(2,2,1);
        barplot(D1.seqNum,D1.MT,'facecolor',{blue},'leg','auto','leglocation','northeast');
        title('Movement time'); xlabel('Sequence'); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); xticklabels(seqlabels);
        subplot(2,2,2);
        [x,~]=lineplot([D2.seqNum D2.day],D2.MT,'errorbars','shade','shadecolor',{blue},'linewidth',lw,'linecolor',{blue},'leg','auto','leglocation','northeast');
        title('Movement time over days'); xlabel('Sequence / Days'); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); xticks(x(2:4:end)); xticklabels(seqlabels);
        % ER
        subplot(2,2,3);
        barplot(D1.seqNum,D1.pER*100,'facecolor',{red},'leg','auto','leglocation','northeast');
        title('Error rate'); xlabel('Sequence'); ylabel('Error rate (%)'); set(gca,'fontsize',fs); xticklabels(seqlabels);
        subplot(2,2,4);
        [x,~]=lineplot([D2.seqNum D2.day],D2.pER*100,'errorbars','shade','shadecolor',{red},'linewidth',lw,'linecolor',{red},'leg','auto','leglocation','northeast');
        title('Error rate over days'); xlabel('Sequence / Days'); ylabel('Error rate (%)'); set(gca,'fontsize',fs); xticks(x(2:4:end)); xticklabels(seqlabels);
        % out
        D.D1=D1; D.D2=D2; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'IPI' % analysis of the inter-press intervals (IPIs) as a function of training and prep time
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % add IPI info
        D.IPI=diff([D.pressTime1,D.pressTime2,D.pressTime3,D.pressTime4,D.pressTime5],1,2);
        D.IPI_1=D.IPI(:,1); D.IPI_2=D.IPI(:,2); D.IPI_3=D.IPI(:,3); D.IPI_4=D.IPI(:,4);
        
%         %%
%         %---------------------------------------------------------------------------------------------------
%         % create summary table for MT
%         T=tapply(D,{'SN','day','train'},...
%             {D.IPI_1,'nanmean','name','IPI1'},...
%             {D.IPI_2,'nanmean','name','IPI2'},...
%             {D.IPI_3,'nanmean','name','IPI3'},...
%             {D.IPI_4,'nanmean','name','IPI4'},...
%             'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.prepTime==2400);
%         for i=1:size(D.IPI,2)
%             T.IPI(:,i)=eval(sprintf('T.IPI%d',i));
%             T=rmfield(T,sprintf('IPI%d',i));
%             T.IPInum(:,i)=ones(size(T.SN,1),1)*i;
%             T.SN(:,i)=T.SN(:,1);
%             T.day(:,i)=T.day(:,1);
%             T.train(:,i)=T.train(:,1);
%         end
%         T.IPI=reshape(T.IPI,size(T.IPI,1)*size(T.IPI,2),1);
%         T.IPInum=reshape(T.IPInum,size(T.IPInum,1)*size(T.IPInum,2),1);
%         T.SN=reshape(T.SN,size(T.SN,1)*size(T.SN,2),1);
%         T.day=reshape(T.day,size(T.day,1)*size(T.day,2),1);
%         T.train=reshape(T.train,size(T.train,1)*size(T.train,2),1);
%         
%         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
%         T = normData(T, {'IPI'}, 'sub');
%         
%         % open figure
%         if nargin>1; figure('Name',sprintf('IPI - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI - group (N=%d)',ns)); end
%         set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
%         
%         % ---------------------------------------------------------------------------------------------------
%         subplot(1,2,1); title('Untrained');
%         plt.line(T.IPInum,T.normIPI, 'split',T.day, 'errorbars','shade', 'style',ipisty, 'leg',dleg, 'leglocation','northeast', 'subset',T.train==0);
%         xlabel('Transition number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
%         
%         subplot(1,2,2); title('Trained');
%         plt.line(T.IPInum,T.normIPI, 'split',T.day, 'errorbars','shade', 'style',ipisty, 'leg',dleg, 'leglocation','northeast', 'subset',T.train==1);
%         xlabel('Transition number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
%         %%
        
        %---------------------------------------------------------------------------------------------------
        % create summary table for MT
        T=tapply(D,{'SN','prepTime','train'},...
            {D.IPI_1,'nanmean','name','IPI1'},...
            {D.IPI_2,'nanmean','name','IPI2'},...
            {D.IPI_3,'nanmean','name','IPI3'},...
            {D.IPI_4,'nanmean','name','IPI4'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.day==4 & D.BN_day>10);
        for i=1:size(D.IPI,2)
            T.IPI(:,i)=eval(sprintf('T.IPI%d',i));
            T=rmfield(T,sprintf('IPI%d',i));
            T.IPInum(:,i)=ones(size(T.SN,1),1)*i;
            T.SN(:,i)=T.SN(:,1);
            T.prepTime(:,i)=T.prepTime(:,1);
            T.train(:,i)=T.train(:,1);
        end
        T.IPI=reshape(T.IPI,size(T.IPI,1)*size(T.IPI,2),1);
        T.IPInum=reshape(T.IPInum,size(T.IPInum,1)*size(T.IPInum,2),1);
        T.SN=reshape(T.SN,size(T.SN,1)*size(T.SN,2),1);
        T.prepTime=reshape(T.prepTime,size(T.prepTime,1)*size(T.prepTime,2),1);
        T.train=reshape(T.train,size(T.train,1)*size(T.train,2),1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'IPI'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.IPInum, T.prepTime, T.normIPI, 'length');
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % ---------------------------------------------------------------------------------------------------
        subplot(2,2,1); title('Untrained');
        plt.line(T.IPInum,T.normIPI, 'split',T.prepTime, 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==0);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        subplot(2,2,2); title('Trained');
        plt.line(T.IPInum,T.normIPI, 'split',T.prepTime, 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==1);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        % stats
        T2 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==1 & (T.prepTime==800 | T.prepTime==1600));
        ttest(T2.normIPI(T2.prepTime==800), T2.normIPI(T2.prepTime==1600), 2, 'paired');
        T3 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==2 & (T.prepTime==400 | T.prepTime==800));
        ttest(T3.normIPI(T3.prepTime==400), T3.normIPI(T3.prepTime==800), 2, 'paired');
        T4 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==4 & (T.prepTime==400 | T.prepTime==2400));
        ttest(T4.normIPI(T4.prepTime==400), T4.normIPI(T4.prepTime==2400), 2, 'paired');
        T5 = tapply(T, {'SN','train'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==4 & T.prepTime==2400);
        ttest(T5.normIPI(T5.train==0), T5.normIPI(T5.train==1), 2, 'paired');
        
        %---------------------------------------------------------------------------------------------------
        % IPIs Comparison
        T2=tapply(D,{'SN','prepTime','train'},...
            {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
            {D.IPI_3,'nanmean','name','IPI3'},...
            {D.IPI_4,'nanmean','name','IPI4'},...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'first2ipi', 'last2ipi', 'IPI3', 'IPI4', 'realPrepTime'}, 'sub');
        
        gap = 2500; % how much space between IPI pair plots (should be > 2400)
        T3.IPI = [T2.normfirst2ipi; T2.normlast2ipi];
        T3.RPT = [T2.normrealPrepTime; T2.normrealPrepTime+gap];
        T3.PT = [T2.prepTime; T2.prepTime+gap];
        T3.train = [T2.train; T2.train];
        T3.pair = [zeros(numel(T2.SN),1); ones(numel(T2.SN),1)] + 1;
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1); title('');
        plt.xy(T3.RPT,T3.IPI,T3.PT, 'split',[T3.pair T3.train], 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); %axis square;
        xticks(unique(T3.PT));
        xticklabels(repmat(unique(T.prepTime)',1,2)); %xticklabels(repmat(unique(T.prepTime)',1,3));
        xlim([min(T3.PT)-100 max(T3.PT)+100]);
        ylim([251 849]); %ylim([300 900]);
        drawline(unique(T3.PT),'dir','vert','linestyle',':','color','k');
        
        %---------------------------------------------------------------------------------------------------
        % stats
        T4=tapply(D,{'SN','prepTime','train'},...
            {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        ttest(T4.first2ipi(T4.prepTime==1600 & T4.train==0), T4.first2ipi(T2.prepTime==1600 & T4.train==1), 2, 'paired');
        ttest(T4.last2ipi(T4.prepTime==1600 & T4.train==0), T4.first2ipi(T2.prepTime==1600 & T4.train==1), 2, 'paired');
        anovaMixed(T4.last2ipi, T4.SN, 'within', [T4.prepTime, T4.train], {'prepTime','train'});
        
        % out
        D.T=T; D.T2=T2; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'force' % data extraction of the force of the presses during sequence execution (peak press force, delta time press-release, pre-press baseline force, peak press velocity, peak release velocity)
        if nargin>1 % load single subj data
            subj=varargin(1);
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj{1}))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % load force data
        b=0; %initialize block number
        for t=1:numel(D.TN)
            if D.BN(t,1)~=b
                b=D.BN(t,1); %update block number within the loop
                MOV=movload(fullfile(pathToData,sprintf('sr2_s%02d_%02d.mov',D.SN(t,1),D.BN(t,1)))); %load MOV file
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
            save(fullfile(pathToAnalyze,sprintf('sr2_force_%s.mat',subj{1})),'-struct', 'D'); %save .mat file
        else % group data
            save(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',numel(subj))),'-struct', 'D'); %save .mat file
        end
        varargout={D}; %return main structure
        
    case 'peak_force' % peak press force
        if nargin>1 % load single subj data
            subj=varargin(1);
            s=str2double(subj{1}(2:end));
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
            D=getrow(D,D.SN==s);
        else % load group data
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
        end
        
        % create summary tables for peak press force
        %         D1=tapply(D,{'SN','BN','day'},...
        %             {D.peak_press_force_T,'nanmedian','name','forceT'},...
        %             {D.peak_press_force_I,'nanmedian','name','forceI'},...
        %             {D.peak_press_force_M,'nanmedian','name','forceM'},...
        %             {D.peak_press_force_R,'nanmedian','name','forceR'},...
        %             {D.peak_press_force_L,'nanmedian','name','forceL'},...
        %             'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','BN','day','train'},...
            {D.peak_press_force_firstPress,'nanmedian','name','forceFirst'},...
            {D.peak_press_force_avgAllFingers,'nanmedian','name','forceAvg'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        D2 = normData(D2, {'forceAvg'}, 'sub');
        
        % open figure
        if nargin>1; figure('Name',sprintf('Peak press force - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Peak press force - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % MT
        %subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        %plt.line([D1.day D1.BN],[D1.forceT,D1.forceI,D1.forceM,D1.forceR,D1.forceL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        %xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); %axis square;
        %subplot(3,1,2); title('Avg all fingers');
        subplot(2,2,1);
        plt.line([D2.day D2.BN],D2.normforceAvg,'split',D2.train,'errorbars','shade','style',trsty_cbs,'leg',trleg,'leglocation','northwest');
        xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
        ylim([1.8 5.2]);
        xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        xt=xticks; drawline(xt([1 12, 13 24, 25 36, 37 48]),'dir','vert','linestyle',':','color','k');
        %xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
        %subplot(3,1,3); title('First press only');
        %plt.line([D2.day D2.BN],D2.forceFirst,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        %xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); %axis square;
        
        % out
        varargout={D}; %return main structure
        
    case 'press_time' % delta time press-release
        if nargin>1 % load single subj data
            subj=varargin(1);
            s=str2double(subj{1}(2:end));
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
            D=getrow(D,D.SN==s);
        else % load group data
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
        end
        
        % create summary tables for peak press force
        D1=tapply(D,{'SN','BN','day'},...
            {D.delta_t_pressRelease_T,'nanmean','name','pressTimeT'},...
            {D.delta_t_pressRelease_I,'nanmean','name','pressTimeI'},...
            {D.delta_t_pressRelease_M,'nanmean','name','pressTimeM'},...
            {D.delta_t_pressRelease_R,'nanmean','name','pressTimeR'},...
            {D.delta_t_pressRelease_L,'nanmean','name','pressTimeL'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','BN','day','train'},...
            {D.delta_t_pressRelease_firstPress,'nanmean','name','pressTimeFirst'},...
            {D.delta_t_pressRelease_avgAllFingers,'nanmean','name','pressTimeAvg'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Delta time press-release - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Delta time press-release - group (N=%d)',ns)); end
        % MT
        subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        plt.line([D1.day D1.BN],[D1.pressTimeT,D1.pressTimeI,D1.pressTimeM,D1.pressTimeR,D1.pressTimeL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Delta time press-release (ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,2); title('Avg all fingers'); sty=style.custom(trc,'markersize',ms,'linewidth',lw);
        plt.line([D2.day D2.BN],D2.pressTimeAvg,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Delta time press-release (ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,3); title('First press only');
        plt.line([D2.day D2.BN],D2.pressTimeFirst,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Delta time press-release (ms)'); set(gca,'fontsize',fs); %axis square;
        
        % out
        varargout={D}; %return main structure
        
    case 'baseline_force' % pre-press baseline force
        if nargin>1 % load single subj data
            subj=varargin(1);
            s=str2double(subj{1}(2:end));
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
            D=getrow(D,D.SN==s);
        else % load group data
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
        end
        
        % create summary tables for peak press force
        D1=tapply(D,{'SN','BN','day'},...
            {D.baseline_force_T,'nanmean','name','baselineT'},...
            {D.baseline_force_I,'nanmean','name','baselineI'},...
            {D.baseline_force_M,'nanmean','name','baselineM'},...
            {D.baseline_force_R,'nanmean','name','baselineR'},...
            {D.baseline_force_L,'nanmean','name','baselineL'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','BN','day','train'},...
            {D.baseline_force_firstPress,'nanmean','name','baselineFirst'},...
            {D.baseline_force_avgAllFingers,'nanmean','name','baselineAvg'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Pre-press baseline force - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Pre-press baseline force - group (N=%d)',ns)); end
        % MT
        subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        plt.line([D1.day D1.BN],[D1.baselineT,D1.baselineI,D1.baselineM,D1.baselineR,D1.baselineL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Baseline force (0-3100 ms; N)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,2); title('Avg all fingers'); sty=style.custom(trc,'markersize',ms,'linewidth',lw);
        plt.line([D2.day D2.BN],D2.baselineAvg,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Baseline force (0-3100 ms; N)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,3); title('First press only');
        plt.line([D2.day D2.BN],D2.baselineFirst,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Baseline force (0-3100 ms; N)'); set(gca,'fontsize',fs); %axis square;
        
        % out
        varargout={D}; %return main structure
        
    case 'peak_velocity' % peak press and release velocity
        if nargin>1 % load single subj data
            subj=varargin(1);
            s=str2double(subj{1}(2:end));
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
            D=getrow(D,D.SN==s);
        else % load group data
            D=load(fullfile(pathToAnalyze,sprintf('sr2_force_group_data_n=%02d.mat',ns)));
        end
        
        % create summary tables for peak press velocity
        D1=tapply(D,{'SN','BN','day'},...
            {D.peak_press_velocity_T,'nanmean','name','pressVelT'},...
            {D.peak_press_velocity_I,'nanmean','name','pressVelI'},...
            {D.peak_press_velocity_M,'nanmean','name','pressVelM'},...
            {D.peak_press_velocity_R,'nanmean','name','pressVelR'},...
            {D.peak_press_velocity_L,'nanmean','name','pressVelL'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','BN','day','train'},...
            {D.peak_press_velocity_firstPress,'nanmean','name','pressVelFirst'},...
            {D.peak_press_velocity_avgAllFingers,'nanmean','name','pressVelAvg'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Peak press velocity - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Peak press velocity - group (N=%d)',ns)); end
        % MT
        subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        plt.line([D1.day D1.BN],[D1.pressVelT,D1.pressVelI,D1.pressVelM,D1.pressVelR,D1.pressVelL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak press velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,2); title('Avg all fingers'); sty=style.custom(trc,'markersize',ms,'linewidth',lw);
        plt.line([D2.day D2.BN],D2.pressVelAvg,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak press velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,3); title('First press only');
        plt.line([D2.day D2.BN],D2.pressVelFirst,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak press velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        
        % create summary tables for peak release velocity
        D1=tapply(D,{'SN','BN','day'},...
            {D.peak_release_velocity_T,'nanmean','name','relVelT'},...
            {D.peak_release_velocity_I,'nanmean','name','relVelI'},...
            {D.peak_release_velocity_M,'nanmean','name','relVelM'},...
            {D.peak_release_velocity_R,'nanmean','name','relVelR'},...
            {D.peak_release_velocity_L,'nanmean','name','relVelL'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        D2=tapply(D,{'SN','BN','day','train'},...
            {D.peak_release_velocity_firstPress,'nanmean','name','relVelFirst'},...
            {D.peak_release_velocity_avgAllFingers,'nanmean','name','relVelAvg'},...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Peak release velocity - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Peak release velocity - group (N=%d)',ns)); end
        % MT
        subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        plt.line([D1.day D1.BN],[D1.relVelT,D1.relVelI,D1.relVelM,D1.relVelR,D1.relVelL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak release velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,2); title('Avg all fingers'); sty=style.custom(trc,'markersize',ms,'linewidth',lw);
        plt.line([D2.day D2.BN],D2.relVelAvg,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak release velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        subplot(3,1,3); title('First press only');
        plt.line([D2.day D2.BN],D2.relVelFirst,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak release velocity (N/ms)'); set(gca,'fontsize',fs); %axis square;
        
        % out
        varargout={D}; %return main structure
        
    case 'press_dur_ET' % analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
        if nargin > 1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr2_%s.mat',subj)) ); %load data for this subject
            D.SN = ones( numel(D.TN), 1) * str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr2_training_all_data.mat') );
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Press duration / overlap - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Press duration / overlap - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %---------------------------------------------------------------------------------------------------
        %D.rel2press_dur(D.rel2press_dur<=0) = NaN; % remove overlap transitions
        % create summary table
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {(sum(D.press2rel_dur,2)./D.MT)*100,'nanmean', 'name','p2r'},... {(D.press2rel_dur,'nanmean', 'name','p2r'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T1.p2r           = reshape(T.p2r, [], 1);
        T1.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T1.BN            = repmat(T.BN, numel(T.p2r(1,:)), 1);
        T1.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T1.day           = repmat(T.day, numel(T.p2r(1,:)), 1);
        T1.P             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T1 = normData(T1, {'p2r'}, 'sub');
        
        subplot(2,2,1); title('');
        plt.bar([T1.day], T1.normp2r, 'split',T1.train, 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        
        %---------------------------------------------------------------------------------------------------
        %D.rel2press_dur(D.rel2press_dur<=0) = NaN; % remove overlap transitions
        % create summary table
        T = tapply(D,{'SN', 'day', 'train'},...
            {D.press2rel_dur./D.MT*100,'nanmean', 'name','p2r'},... {D.press2rel_dur,'nanmean', 'name','p2r'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T2.p2r           = reshape(T.p2r, [], 1);
        T2.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T2.day           = repmat(T.day, numel(T.p2r(1,:)), 1);
        T2.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T2.P             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'p2r'}, 'sub');
        
        subplot(2,2,2); title('');
        %plt.line(T2.P, T2.normp2r, 'plotfcn','nanmean', 'split',T2.day, 'errorbars','shade', 'style',ipisty, 'leg',dleg, 'leglocation','northeast');
        plt.line(T2.P, T2.normp2r, 'plotfcn','nanmean', 'split',T2.day, 'errorbars','shade', 'style',ptUNsty, 'leg',dleg, 'leglocation','northeast', 'subset',T2.train==0);
        hold on;
        plt.line(T2.P, T2.normp2r, 'plotfcn','nanmean', 'split',T2.day, 'errorbars','shade', 'style',ptTRsty, 'leg',dleg, 'leglocation','northeast', 'subset',T2.train==1);
        xlabel('Press number'); ylabel('Press duration (% of ET)'); set(gca,'fontsize',fs); axis square;
        
        %---------------------------------------------------------------------------------------------------
        %D.rel2press_dur(D.rel2press_dur<=0) = NaN; % remove overlap transitions
        % create summary table
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {(sum(D.rel2press_dur,2)./D.MT)*100,'nanmean', 'name','r2p'},... {D.rel2press_dur,'nanmean', 'name','r2p'},...
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
        %[x,~,~]=plt.line([T3.day T3.BN], T3.normr2p, 'plotfcn','nanmean', 'split',T3.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg, 'leglocation','east');
        %plt.bar([T3.day T3.BN], T3.normr2p, 'split',T3.train, 'style',trsty_cbs, 'leg',trleg, 'leglocation','east');
        plt.box([[T1.day;T3.day], [T1.BN;T3.BN]], [-T1.normp2r; T3.normr2p], 'split',[T1.train; T3.train], 'style',trsty_cbs, 'leg','skip', 'plotall',2, 'whiskerlength',0);
        %xlabel('Block number'); ylabel('Gap duration (% of ET)'); set(gca,'fontsize',fs); axis square;
        xlabel('Block number'); ylabel('Press duration   (% of ET)   Gap duration'); set(gca,'fontsize',fs); axis square;
        x=xticks; xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        hold on; yticklabels([150 100 50 0 50 100 150]);
        
        %---------------------------------------------------------------------------------------------------
        %D.rel2press_dur(D.rel2press_dur<=0) = NaN; % remove overlap transitions
        % create summary table
        T = tapply(D,{'SN', 'day', 'train'},...
            {(D.rel2press_dur./D.MT)*100,'nanmean', 'name','r2p'},... {D.rel2press_dur,'nanmean', 'name','r2p'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T4.r2p           = reshape(T.r2p, [], 1);
        T4.SN            = repmat(T.SN, numel(T.r2p(1,:)), 1);
        T4.day           = repmat(T.day, numel(T.r2p(1,:)), 1);
        T4.train         = repmat(T.train, numel(T.r2p(1,:)), 1);
        T4.T             = reshape(repmat(1:numel(T.r2p(1,:)), numel(T.r2p(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T4 = normData(T4, {'r2p'}, 'sub');
        
        subplot(2,2,4); title('');
        %plt.line(T4.T, T4.normr2p, 'plotfcn','nanmean', 'split',T4.day, 'errorbars','shade', 'style',ipisty, 'leg',dleg, 'leglocation','northeast');
        plt.line(T4.T, T4.normr2p, 'plotfcn','nanmean', 'split',T4.day, 'errorbars','shade', 'style',ptUNsty, 'leg',dleg, 'leglocation','northeast', 'subset',T4.train==0);
        hold on;
        plt.line(T4.T, T4.normr2p, 'plotfcn','nanmean', 'split',T4.day, 'errorbars','shade', 'style',ptTRsty, 'leg',dleg, 'leglocation','northeast', 'subset',T4.train==1);
        xlabel('Transition number'); ylabel('Gap duration (% of ET)'); set(gca,'fontsize',fs); axis square;
        
        hold on;
        subplot(2,2,1);
        ylim([-139 139]);
        drawline(x([1 12, 13 24, 25 36, 37 48]),'dir','vert','linestyle',':','color','k');
        drawline(0, 'dir','horz', 'linestyle','-', 'color','k');
        drawline(100, 'dir','horz', 'linestyle',':', 'color','k');
        drawline(-100, 'dir','horz', 'linestyle',':', 'color','k');
        subplot(2,2,2);
        subplot(2,2,4);
        drawline(0, 'dir','horz', 'linestyle',':', 'color','k');
        hold off;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'press_dur_ms' % analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
        if nargin > 1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr2_%s.mat',subj)) ); %load data for this subject
            D.SN = ones( numel(D.TN), 1) * str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr2_training_all_data.mat') );
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Press duration / overlap - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Press duration / overlap - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %---------------------------------------------------------------------------------------------------
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
        
        subplot(2,2,2); title(''); hold on;
        plt.line([T1.day T1.BN], T1.normp2r, 'plotfcn','nanmean', 'split',T1.train, 'errorbars','shade', 'style',dursty, 'leg',trleg, 'leglocation','southeast');
        xlabel('Block number'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        xticklabels({'3','','','','','','','','','','','',... '14',...
            '17','','','','','','','','','','','',... '28',...
            '31','','','','','','','','','','','',... '42',...
            '45','','','','','','','','','','','56'});
        
        %---------------------------------------------------------------------------------------------------
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
        
    case 'press_dur_PT' % analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
        if nargin > 1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr2_%s.mat',subj)) ); %load data for this subject
            D.SN = ones( numel(D.TN), 1) * str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr2_training_all_data.mat') );
        end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Press duration / overlap - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Press duration / overlap - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %---------------------------------------------------------------------------------------------------
        % create summary table
        T = tapply(D,{'SN', 'train', 'prepTime'},...
            {sum(D.press2rel_dur,2),'nanmean', 'name','p2r'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0 & D.day==4);
        T1.p2r           = reshape(T.p2r, [], 1);
        T1.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T1.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T1.prepTime      = repmat(T.prepTime, numel(T.p2r(1,:)), 1);
        T1.P             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T1 = normData(T1, {'p2r'}, 'sub');
        
        subplot(2,2,1); title('Day 4');
        plt.line([T1.prepTime], T1.normp2r, 'plotfcn','nanmean', 'split',T1.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        %---------------------------------------------------------------------------------------------------
        % create summary table
        T = tapply(D,{'SN', 'prepTime', 'train'},...
            {D.press2rel_dur,'nanmean', 'name','p2r'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0 & D.day==4);
        T2.p2r           = reshape(T.p2r, [], 1);
        T2.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T2.prepTime      = repmat(T.prepTime, numel(T.p2r(1,:)), 1);
        T2.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T2.P             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'p2r'}, 'sub');
        
        subplot(2,2,2); title('');
        plt.line(T2.P, T2.normp2r, 'plotfcn','nanmean', 'split',T2.prepTime, 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T2.train==0);
        hold on;
        plt.line(T2.P, T2.normp2r, 'plotfcn','nanmean', 'split',T2.prepTime, 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T2.train==1);
        xlabel('Press number'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        %---------------------------------------------------------------------------------------------------
        % create summary table
        T = tapply(D,{'SN', 'train', 'prepTime'},...
            {sum(D.rel2press_dur,2),'nanmean', 'name','r2p'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0 & D.day==4);
        T3.r2p           = reshape(T.r2p, [], 1);
        T3.SN            = repmat(T.SN, numel(T.r2p(1,:)), 1);
        T3.train         = repmat(T.train, numel(T.r2p(1,:)), 1);
        T3.prepTime      = repmat(T.prepTime, numel(T.r2p(1,:)), 1);
        T3.T             = reshape(repmat(1:numel(T.r2p(1,:)), numel(T.r2p(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T3 = normData(T3, {'r2p'}, 'sub');
        
        subplot(2,2,3); title('Day 4');
        plt.line(T3.prepTime, T3.normr2p, 'plotfcn','nanmean', 'split',T3.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time'); ylabel('Gap duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        %---------------------------------------------------------------------------------------------------
        % create summary table
        T = tapply(D,{'SN', 'prepTime', 'train'},...
            {D.rel2press_dur,'nanmean', 'name','r2p'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0 & D.day==4);
        T4.r2p           = reshape(T.r2p, [], 1);
        T4.SN            = repmat(T.SN, numel(T.r2p(1,:)), 1);
        T4.prepTime      = repmat(T.prepTime, numel(T.r2p(1,:)), 1);
        T4.train         = repmat(T.train, numel(T.r2p(1,:)), 1);
        T4.T             = reshape(repmat(1:numel(T.r2p(1,:)), numel(T.r2p(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T4 = normData(T4, {'r2p'}, 'sub');
        
        subplot(2,2,4); title('');
        plt.line(T4.T, T4.normr2p, 'plotfcn','nanmean', 'split',T4.prepTime, 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T4.train==0);
        hold on;
        plt.line(T4.T, T4.normr2p, 'plotfcn','nanmean', 'split',T4.prepTime, 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T4.train==1);
        xlabel('Transition number'); ylabel('Gap duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        hold on;
        subplot(2,2,1);
        subplot(2,2,2);
        subplot(2,2,3);
        drawline(0, 'dir','horz', 'linestyle',':', 'color','k');
        subplot(2,2,4);
        drawline(0, 'dir','horz', 'linestyle',':', 'color','k');
        hold off;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 'rtt_rpt' % analysis of reaction time task (rtt), with real preparation time (rpt)
        
        if nargin>1 % load single subj data
            subj=varargin{1};
            D=load(fullfile(pathToData,sprintf('sr2_%s.mat',subj))); %load data for this subject
            D.SN=ones(numel(D.TN),1)*str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D=load(fullfile(pathToAnalyze,'sr2_training_all_data.mat'));
        end
        
        % ----------------------------------------------------------------------------------------------------------------------------
        % create summary tables for ACC
        T = tapply(D, {'SN','day','prepTime'}, ...
            {(1-D.isError)*100,'nanmean','name','ACC'}, ... {(1-D.pressError)*100,'nanmean','name','ACC'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','prepTimeReal'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==1 & D.timingError==0);
        
        % normalize ACC data to remove between-subject variability
        T = normData(T, {'ACC', 'prepTimeReal'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.prepTime, T.day, T.normACC, 'length');
        
        % open figure
        if nargin>1; figure('Name',sprintf('RTT prepTimeReal - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('RTT prepTimeReal - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1);
        [~,~]=plt.xy(T.normprepTimeReal,T.normACC,T.prepTime, 'errorbars','plusminus_wocap', 'split',T.day, 'style',ipisty, 'leg',dleg, 'leglocation','southeast');
        xlim([175 675]); ylim([16 104]);
        xlabel('True reaction time (ms)'); ylabel('Finger selection accuracy (%)');
        axis square; set(gca,'fontsize',fs);
        drawline(20,'dir','horz','linestyle','-','color',black);
        hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        % ----------------------------------------------------------------------------------------------------------------------------
        % log fit
        rtt_x_day=zeros(numel(unique(T.day)),1);
        rtt_y=0.80; % the accuracy value (from 0 to 1) that you want to use to determine subj RT with that accuracy (from logistic fit)
        % fit logistic function
        R.x=[]; R.y_hat=[]; R.day=[]; R.prepTime=[];
        Z = getrow(D, D.dummy==0 & D.rtt==1 & D.timingError==0);
        for iDay = 1:numel(unique(Z.day))
            xpt = Z.prepTime(Z.day==iDay) + Z.RT(Z.day==iDay);
            y = ~Z.isError(Z.day==iDay);
            theta_zero = [0.01, 10]';
            [theta_hat] = fitlog(xpt, y, theta_zero);
            a = theta_hat(1); b = theta_hat(2);
            c = .20;
            rtt_x_day(iDay) = (b - log((1-c)/(rtt_y-c) - 1)) / a;
            x = (150:50:700)';
            y_hat = modlog(theta_hat, x);
            R.x = [R.x; x];
            R.y_hat = [R.y_hat; y_hat];
            R.day = [R.day; ones(numel(x),1)*iDay];
            R.prepTime = [R.prepTime; Z.prepTime(Z.day==iDay)];
        end
        
        % out
        D.T=T; D.rtt_y=rtt_y; D.rtt_x_day=rtt_x_day; D.R=R; %incorporate the sub-structures as fields of main structure
        varargout={D,rtt_x_day}; %return main structure
        
    case 'rtt_rpt_logfit' % analysis of reaction time task (rtt), with real preparation time (rpt), per subject (subj), including logfit
        
        if nargin>1 % load single subj data
            subj=varargin;
            ns=numel(subj);
            subvec=zeros(1,ns); for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end
        else % load group data
            subj={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13',	's15','s16','s17','s18','s19','s20'};
            ns=numel(subj);
            subvec=zeros(1,ns); for i=1:ns; subvec(1,i)=str2double(subj{i}(2:3)); end
        end
        
        %---------------------------------------------------------------------------------------------------
        Z.x = []; Z.y = []; Z.day = []; Z.pt = [];
        Z.SN = []; rtt_x_day = zeros(ns, 4);
        for s=1:ns
            sn=subj{s};
            S=sr2_training_analyze('rtt_rpt',sn); close all;
            rtt_x_day(s,:)=real(S.rtt_x_day);
            % store info about model fits per subject
            Z.x      = [Z.x; S.R.x];
            Z.y      = [Z.y; S.R.y_hat];
            Z.day    = [Z.day; S.R.day];
            Z.pt     = [Z.pt; S.R.prepTime];
            Z.SN     = [Z.SN; ones(numel(S.R.x),1) * subvec(s)];
        end
        Z.rtt_x_all = rtt_x_day;
        Z.rtt_x = mean(rtt_x_day, 1);
        Z.rtt_y = S.rtt_y;
        save( fullfile(pathToAnalyze, 'sr2_model_fits_subj_rpt.mat'), '-struct', 'Z'); %save all_data.mat file
        
        %---------------------------------------------------------------------------------------------------
        % produce old plot
        sr2_training_analyze('rtt_rpt'); hold on;
        
        subplot(2,2,3); title('Logistic fit'); hold on;
        plt.line(Z.x,Z.y*100, 'plotfcn','nanmean', 'split',Z.day, 'style',rttModelsty, 'leg',dleg, 'leglocation','southeast');
        xticks(200:50:650); xticklabels({'200','','300','','400','','500','','600',''});
        xlim([175 675]); ylim([16 104]);
        xlabel('Preparation time (ms)'); ylabel('Finger selection accuracy (%)');
        axis square; set(gca,'fontsize',fs);
        hold on;
        drawline(80,'dir','horz','linestyle','-','color',black);
        drawline(400,'dir','vert','linestyle','-','color',black);
        hold off;
        
        %---------------------------------------------------------------------------------------------------
        % stats
        % gain in selection accuracy
        ST = tapply(Z,{'SN', 'day'},...
            {Z.y*100, 'nanmean', 'name','y', 'subset',Z.x==400});
        ttest(ST.y(ST.day==1), ST.y(ST.day==4), 2, 'paired');
        
        % gain in planning speed
        ttest(Z.rtt_x_all(:, 1), Z.rtt_x_all(:, 4), 2, 'paired');
        
        % out
        varargout={Z}; %return main structure
        
    otherwise
        error('no such case!')
end
end

function [theta_hat]=fitlog(x,y,theta_zero)
fcn=@(theta)bernoulli_loss(y,modlog(theta,x));
theta_hat=fminsearch(fcn,theta_zero);
end

function [y_hat]=modlog(theta,x)
a=theta(1,1); b=theta(2,1); %c=theta(3,1);
c=.20; %c=0;
y_hat=1./(1 + exp(-a*x + b))*(1-c) + c; %y_hat=1./(1 + exp(-a*x + b));
end

function [loss]=bernoulli_loss(y,y_hat)
n=numel(y);
loss=zeros(n,1);
for i=1:n; loss(i,1)=(y(i)*log(y_hat(i)) + (1-y(i))*log(1-y_hat(i))); end
loss=-(nansum(loss));
end
