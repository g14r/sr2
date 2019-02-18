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
%                       [D]=sr2_training_analyze('rtt_logfit_subj');                    %single subject results of logistic fit estimation (given a certain ACC/y, the corresponding RT/x)
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
%                       [D]=sr2_training_analyze('press_dur');                          %group          analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
%                       [D]=sr2_training_analyze('press_dur','s01');                    %single subject analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
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
darkgray = [80,80,80]/255;
gray3 = [50,50,50]/255;
gray2 = [100,100,100]/255; 
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
ptc_un_cbs={[229	158	102]/255;[221	126	50]/255;[213	94	0]/255;[149	65	0]/255};
ptc_tr_cbs={[153	198	224]/255;[76	156	201]/255;[0	102	160]/255;[0	68	106]/255};
% ptc_un={[254,204,92]/255;[253,141,60]/255;[240,59,32]/255;[189,0,38]/255};
% ptc_tr={[194,230,153]/255;[120,198,121]/255;[49,163,84]/255;[0,104,55]/255};
% c_rtt_data={[200,200,200]/255,[150,150,150]/255,[100,100,100]/255,[0,0,0]/255};
% c_rtt_model={[158,202,225]/255,[107,174,214]/255,[33,113,181]/255,[8,69,148]/255};
%c_rtt_data_cbs={[102	196	171]/255;[0	158	115]/255;[	0	107	78]/255;[0	82	59]/255};
%c_rtt_model_cbs={[219	161	193]/255;[204	121	167]/255;[163	96	133]/255;[102	60	83]/255};
fingc={green,yellow,red,purple,blue};
trc={lightgray,black};

% plot defaults
fs=20; %default font size for all figures
lw=2;%3; %4; %3; %default line width for all figures
ms=8;%12; %10; %default marker size for all figures

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,gray,lightgray,green,lightgreen,black,silver,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
trsty_cbs=style.custom({cbs_red,cbs_blue},'markersize',ms,'linewidth',lw,'markersize',ms,'errorwidth',lw);%,'markertype','none');
unday1sty_cbs=style.custom({cbs_red,cbs_blue},'markersize',ms,'linewidth',lw,'linestyle',':','markersize',ms,'errorwidth',lw);%,'markertype','none');
ptUNsty=style.custom(ptc_un_cbs,'markersize',ms,'linewidth',lw);
ptTRsty=style.custom(ptc_tr_cbs,'markersize',ms,'linewidth',lw);
%rttDatasty=style.custom(c_rtt_data_cbs,'markersize',ms,'linewidth',lw);
%rttModelsty=style.custom(c_rtt_model_cbs,'markersize',ms,'linewidth',lw);
%darkgraysty=style.custom(black,'markersize',ms,'linewidth',lw,'errorwidth',lw);
ipisty=style.custom({lightgray,gray,darkgray,black},'markersize',ms,'linewidth',lw,'errorwidth',lw);
presssty=style.custom({lightgray,gray,gray2,gray3,black},'markersize',ms,'linewidth',lw,'markersize',ms,'errorwidth',lw);
trnssty=style.custom({gray,gray2,gray3,black},'markersize',ms,'linewidth',lw,'markersize',ms,'errorwidth',lw);

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
        
        rtt_x_day=zeros(numel(unique(T.day)),1);
        rtt_y=0.80; % the accuracy value (from 0 to 1) that you want to use to determine subj RT with that accuracy (from logistic fit)
        % fit logistic function
        R.x=[]; R.y_hat=[]; R.day=[]; y=y/100;
        for iDay=1:numel(unique(T.day))
            %theta_zero=[0.1,400,.20]';
            %theta_zero=[0.01, 5]';
            theta_zero=[0.01, 10]';
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
        %plt.line(R.x,R.y_hat*100,'split',R.day,'style',rttModelsty,'leg',dleg,'leglocation','east');
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
        
    case 'rtt_logfit_subj' % single subject results of logistic fit estimation (given a certain ACC/y, the corresponding RT/x)
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
        
        % produce old plot
        sr2_training_analyze('rtt'); hold on;
        
        % normalize ACC data to remove between-subject variability (i.e. plot within-subject standard error)
        D.SN = D.model_fits_SN;
        D = normData(D, {'model_fits_y'}, 'sub');
        
        subplot(2,2,3);
        %[~,~] = plt.line(D.model_fits_x, D.model_fits_y*100, 'errorbars','shade', 'split',D.model_fits_day, 'style',rttModelsty, 'leg',dleg, 'leglocation','east');
        [~,~] = plt.line(D.model_fits_x, D.model_fits_y*100, 'errorbars','shade', 'split',D.model_fits_day, 'style',ipisty, 'leg',dleg, 'leglocation','east');
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
        %         drawline(D.rtt_x(1), 'dir','vert', 'linestyle','-', 'color',red);
        %         drawline(D.rtt_x(4), 'dir','vert', 'linestyle','-', 'color',red);
        % gain in selection accuracy
        %drawline(400, 'dir','vert', 'linestyle',':', 'color',black);
        %         drawline(y(1,5), 'dir','horz', 'linestyle','-', 'color',red);
        %         drawline(y(4,5), 'dir','horz', 'linestyle','-', 'color',red);
        hold off;
        
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
        
        % create summary tables for MT
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {D.MT,'mean', 'name','MT'},...
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
        %         xticklabels({'3','','','','','Day 1','','','','','','',... '14',...
        %             '17','','','','','Day 2','','','','','','',... '28',...
        %             '31','','','','','Day 3','','','','','','',... '42',...
        %             '45','','','','','Day 4','','','','','','',...'56'
        %             });
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
        
        % stats
        %         ttest(T.ACC(T.day==4 & T.train==0), T.ACC(T.day==4 & T.train==1), 2, 'paired');
        %         [T.t, T.p] = ttest(T.ACC(T.day==4 & T.train==0), T.ACC(T.day==4 & T.train==1), 2, 'paired');
        
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
        
        % stats
        ttest(T.pACC(T.day==4 & T.train==0), T.pACC(T.day==4 & T.train==1), 2, 'paired');
        [T.t, T.p] = ttest(T.pACC(T.day==4 & T.train==0), T.pACC(T.day==4 & T.train==1), 2, 'paired');
        
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
        
        % stats
        %         ttest(T.tACC(T.day==4 & T.train==0), T.tACC(T.day==4 & T.train==1), 2, 'paired');
        %         [T.t, T.p] = ttest(T.tACC(T.day==4 & T.train==0), T.tACC(T.day==4 & T.train==1), 2, 'paired');
        
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
        
        % create summary table for MT
        T = tapply(D, {'SN','prepTime','day'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.SN, T.prepTime, T.MT, 'mean');
        
        % plot data
        subplot(2,2,1); title('');
        %plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'errorbars','plusminus_wocap', 'style',darkgraysty);
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'errorbars','plusminus_wocap', 'split',T.day, 'style',ipisty, 'leg',dleg);
        xlabel('Preparation time (ms)'); ylabel('ET (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([1000 1700]);
        hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        % stats
        anovaMixed(T.MT, T.SN, 'within', [T.prepTime, T.day], {'prepTime','day'});
        T2 = tapply(D, {'SN','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.isError==0 & D.dummy==0 & D.rtt==0 & D.BN_day>10);
        anovaMixed(T2.MT, T2.SN, 'within', [T2.prepTime], {'prepTime'});
        %         ttest(T.MT(T.prepTime==400), T.MT(T.prepTime==800), 2, 'paired');
        %         ttest(T.MT(T.prepTime==800), T.MT(T.prepTime==1600), 2, 'paired');
        %         ttest(T.MT(T.prepTime==1600), T.MT(T.prepTime==2400), 2, 'paired');
        
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
        
        % create summary table for MT (day 4)
        T = tapply(D, {'SN','train','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'MT'}, 'sub');
        
        % compute the percentage of MT at prep time 2400 to show the interaction
        %h=figure; [~,y]=plt.xy(T.realPrepTime,T.normMT,T.prepTime,'split',T.train); close(h);
        h=figure; [~,y]=plt.xy(T.realPrepTime,(1./T.normMT),T.prepTime,'split',T.train); close(h);
        T.pctMT = zeros(numel(T.normMT),1);
        tridx = T.train==1;
        unidx = T.train==0;
        %         T.pctMT(unidx) = 100 + (100 - (T.normMT(unidx) * 100) ./ y(1, numel(unique(T.prepTime))) );
        %         T.pctMT(tridx) = 100 + (100 - (T.normMT(tridx) * 100) ./ y(2, numel(unique(T.prepTime))) );
        T.pctMT(unidx) = (1 ./ (T.normMT(unidx)))  ./  max(y(1,:)) * 100;
        T.pctMT(tridx) = (1 ./ (T.normMT(tridx)))  ./  max(y(2,:)) * 100;
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.train, T.prepTime, T.normMT, 'length'); pivottable(T.train, T.MTdiffPTcat, T.normMTdiffPT, 'length');
        
        % day 1
        %         T2 = tapply(D, {'SN','prepTime'}, ...
        %             {D.MT,'nanmean', 'name','MT'}, ...
        %             {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
        %             'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==1 & D.BN_day>10 & D.train==0 );
        %
        %         % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %         T2 = normData(T2, {'MT'}, 'sub');
        %
        %         % compute the percentage of MT at prep time 2400 to show the interaction
        %         h=figure; [~,y]=plt.xy(T2.realPrepTime,(1./T2.normMT),T2.prepTime); close(h);
        %         T2.pctMT = zeros(numel(T2.normMT),1);
        %         T2.pctMT = (1 ./ (T2.normMT))  ./  max(y(1,:)) * 100;
        
        T2 = tapply(D, {'SN','train','prepTime'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'MT'}, 'sub');
        
        % compute the percentage of MT at prep time 2400 to show the interaction
        h=figure; [~,y]=plt.xy(T2.realPrepTime,(1./T2.normMT),T2.prepTime,'split',T2.train); close(h);
        T2.pctMT = zeros(numel(T2.normMT),1);
        tridx = T2.train==1;
        unidx = T2.train==0;
        T2.pctMT(unidx) = (1 ./ (T2.normMT(unidx)))  ./  max(y(1,:)) * 100;
        T2.pctMT(tridx) = (1 ./ (T2.normMT(tridx)))  ./  max(y(2,:)) * 100;
        
        % plot data
        subplot(2,2,1); title('Day 4 (end of training)');
        plt.xy(T2.realPrepTime,T2.normMT,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs, 'subset',T2.train==0);
        %plt.xy(T2.realPrepTime,T2.normMT,T2.prepTime, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs);
        hold on;
        plt.xy(T.realPrepTime,T.normMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Execution time (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([850 1850]);
        hold on; drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k'); hold off;
        
        subplot(2,2,2); title('Interaction');
        plt.xy(T.realPrepTime,T.pctMT,T.prepTime, 'split',T.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','east');
        hold on;
        plt.xy(T2.realPrepTime,T2.pctMT,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs, 'subset',T2.train==0);
        %plt.xy(T2.realPrepTime,T2.pctMT,T2.prepTime, 'errorbars','plusminus_wocap', 'style',unday1sty_cbs);
        hold on;
        xlabel('Preparation time (ms)'); ylabel('% of ET at 2400'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime(tridx))); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([75 105]);
        hold on; drawline(unique(T.prepTime(tridx)),'dir','vert','linestyle',':','color','k'); drawline(100,'dir','horz','linestyle','--','color','k'); hold off;
        
        % stats
        T.res = anovaMixed(T.MT,T.SN,'within', [T.prepTime,T.train], {'prepTime','train'});
        %ttest(T.MT(T.prepTime==2400 & T.train==0), T.MT(T.prepTime==2400 & T.train==1), 2, 'paired');
        %[T.t, T.p] = ttest(T.MT(T.prepTime==2400 & T.train==0), T.MT(T.prepTime==2400 & T.train==1), 2, 'paired');
        
        T3 = tapply(D, {'SN','prepTime','day'}, ...
            {D.MT,'nanmean', 'name','MT'}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.train==0 & ismember(D.day,[1,4]));
        T3.res = anovaMixed(T3.MT,T3.SN,'within', [T3.prepTime,T3.day], {'prepTime','day'});
        
        % create summary tables for press ACC
        T = tapply(D,{'SN','prepTime','train'},...
            {(1-D.pressError)*100,'nanmean','name','pACC', 'subset',D.timingError==0}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.day==4 & D.BN_day>10);
        
        % normalize pACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'pACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.pACC, 'length');
        
        subplot(2,2,3); hold on; title('Learning - press ACC');
        plt.bar(T.prepTime,T.normpACC, 'split',T.train, 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Training day'); ylabel('Press accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); axis square;
        
        % stats
        %ttest(T.pACC(T.train==0), T.pACC(T.train==1), 2, 'paired');
        %[T.t, T.p] = ttest(T.pACC(T.train==0), T.pACC(T.train==1), 2, 'paired');
        
        % create summary tables for timing ACC
        T = tapply(D,{'SN','prepTime','train'},...
            {(1-D.timingError)*100,'nanmean','name','tACC'}, ...
            'subset',D.dummy==0 & D.rtt==0 & D.day==4 & D.BN_day>10);
        
        % normalize tACC data to remove between-subject variability (i.e. plot within-subject standard error)
        T = normData(T, {'tACC'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T.day, T.BN, T.tACC, 'length');
        
        subplot(2,2,4); hold on; title('Learning - timing ACC');
        plt.bar(T.prepTime,T.normtACC, 'split',T.train, 'style',trsty_cbs, 'leg',{'Untrained','Trained'}, 'leglocation','northeast');
        xlabel('Training day'); ylabel('Timing accuracy (%)'); set(gca,'fontsize',fs); ylim([0 100]); axis square;
        
        % stats
        %ttest(T.tACC(T.train==0), T.tACC(T.train==1), 2, 'paired');
        %[T.t, T.p] = ttest(T.tACC(T.train==0), T.tACC(T.train==1), 2, 'paired');
        
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
        
        subplot(2,2,1); title('Untrained');
        plt.line(T.IPInum,T.normIPI, 'split',T.prepTime, 'errorbars','shade', 'style',ptUNsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==0);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        subplot(2,2,2); title('Trained');
        plt.line(T.IPInum,T.normIPI, 'split',T.prepTime, 'errorbars','shade', 'style',ptTRsty, 'leg',ptleg, 'leglocation','northeast', 'subset',T.train==1);
        xlabel('Interval number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; ylim([50 500]);
        
        % stats
        %         T.anova    = anovaMixed(T.IPI,T.SN,'within', [T.IPInum,T.prepTime], {'IPInum','prepTime'});
        %         T.anova_un = anovaMixed(T.IPI,T.SN,'within', [T.IPInum,T.prepTime], {'IPInum','prepTime'}, 'subset',T.train==0);
        %         T.anova_tr = anovaMixed(T.IPI,T.SN,'within', [T.IPInum,T.prepTime], {'IPInum','prepTime'}, 'subset',T.train==1);
        %
        %         ttest(T.normIPI(T.IPInum==1 & T.train==0 & T.prepTime==800), T.normIPI(T.IPInum==1 & T.train==0 & T.prepTime==1600), 2, 'paired');
        %         ttest(T.normIPI(T.IPInum==1 & T.train==1 & T.prepTime==800), T.normIPI(T.IPInum==1 & T.train==1 & T.prepTime==1600), 2, 'paired');
        %
        %ttest(T.normIPI(T.IPInum==2 & T.train==0 & T.prepTime==400), T.normIPI(T.IPInum==2 & T.train==0 & T.prepTime==800), 2, 'paired');
        %ttest(T.normIPI(T.IPInum==2 & T.train==1 & T.prepTime==400), T.normIPI(T.IPInum==2 & T.train==1 & T.prepTime==800), 2, 'paired');
        %
        %         ttest(T.normIPI(T.IPInum==4 & T.train==0 & T.prepTime==400), T.normIPI(T.IPInum==4 & T.train==0 & T.prepTime==2400), 2, 'paired');
        %ttest(T.normIPI(T.IPInum==4 & T.train==1 & T.prepTime==400), T.normIPI(T.IPInum==4 & T.train==1 & T.prepTime==2400), 2, 'paired');
        
        %         T2 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==1 & (T.prepTime==800 | T.prepTime==1600));
        %         ttest(T2.normIPI(T2.prepTime==800), T2.normIPI(T2.prepTime==1600), 2, 'paired');
        %         T3 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==2 & (T.prepTime==400 | T.prepTime==800));
        %         ttest(T3.normIPI(T3.prepTime==400), T3.normIPI(T3.prepTime==800), 2, 'paired');
        %         T4 = tapply(T, {'SN','prepTime'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==4 & (T.prepTime==400 | T.prepTime==2400));
        %         ttest(T4.normIPI(T4.prepTime==400), T4.normIPI(T4.prepTime==2400), 2, 'paired');
        T5 = tapply(T, {'SN','train'}, {T.normIPI,'nanmean','name','normIPI'}, 'subset',T.IPInum==4 & T.prepTime==2400);
        ttest(T5.normIPI(T5.train==0), T5.normIPI(T5.train==1), 2, 'paired');
        
        % IPIs Comparison
        T2=tapply(D,{'SN','prepTime','train'},...
            {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
            {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
            {D.prepTime + D.RT,'nanmean', 'name','realPrepTime'}, ... % compute actual preparation time
            'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'first2ipi','last2ipi'}, 'sub');
        
        % make sure that you have one value per subject for each condition
        % pivottable(T2.train, T2.prepTime, T2.normfirst2ipi, 'length'); % pivottable(T2.train, T2.prepTime, T2.normlast2ipi, 'length');
        
        % open figure
        if nargin>1; figure('Name',sprintf('IPI comparison - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('IPI comparison - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,3); title('First 2 IPIs');
        plt.xy(T2.realPrepTime,T2.normfirst2ipi,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([300 900]);
        
        subplot(2,2,4); title('Last 2 IPIs');
        plt.xy(T2.realPrepTime,T2.normlast2ipi,T2.prepTime, 'split',T2.train, 'errorbars','plusminus_wocap', 'style',trsty_cbs, 'leg',trleg, 'leglocation','northeast');
        xlabel('Preparation time (ms)'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); axis square; xticks(unique(T.prepTime)); xlim([min(T.prepTime)-100 max(T.prepTime)+100]); ylim([300 900]);
        
        hold on;
        subplot(2,2,3); title('First 2 IPIs');
        drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k');
        subplot(2,2,4); title('Last 2 IPIs');
        drawline(unique(T.prepTime),'dir','vert','linestyle',':','color','k');
        hold off;
        
        % stats
        
        T3.SN = [T2.SN; T2.SN];
        T3.IPI = [T2.first2ipi; T2.last2ipi];
        T3.IPIpair = [ones(numel(T2.first2ipi),1); ones(numel(T2.last2ipi),1)*2];
        T3.train = [T2.train; T2.train];
        T3.prepTime = [T2.prepTime; T2.prepTime];
        T3 = tapply(T3, {'SN','train','IPIpair'}, {T3.IPI,'nanmean', 'name','IPI'}, 'subset',T3.prepTime==2400);
        T.anova = anovaMixed(T3.IPI, T3.SN, 'within', [T3.IPIpair, T3.train], {'IPIpair','train'});
        
        T.anova    = anovaMixed(T2.last2ipi, T2.SN, 'within', [T2.prepTime, T2.train], {'prepTime','train'});
        %         T3=tapply(D,{'SN','train'},...
        %             {(D.IPI_1+D.IPI_2),'nanmean','name','first2ipi'},...
        %             'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10 & D.prepTime==2400);
        %         ttest(T3.first2ipi(T3.train==0), T3.first2ipi(T3.train==1), 2, 'paired');
        %
        %         T5=tapply(D,{'SN','train'},...
        %             {(D.IPI_3+D.IPI_4),'nanmean','name','last2ipi'},...
        %             'subset',D.dummy==0 & D.rtt==0 & D.isError==0 & D.day==4 & D.BN_day>10 & D.prepTime==2400);
        %         ttest(T5.last2ipi(T5.train==0), T5.last2ipi(T5.train==1), 2, 'paired');
        
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
        % MT
        %subplot(3,1,1); title('Five fingers'); sty=style.custom(fingc,'markersize',ms,'linewidth',lw);
        %plt.line([D1.day D1.BN],[D1.forceT,D1.forceI,D1.forceM,D1.forceR,D1.forceL],'errorbars','shade','style',sty,'leg',fingleg,'leglocation','northeast');
        %xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); %axis square;
        %subplot(3,1,2); title('Avg all fingers');
        sty=style.custom(trc,'markersize',ms,'linewidth',lw);
        plt.line([D2.day D2.BN],D2.normforceAvg,'split',D2.train,'errorbars','shade','style',sty,'leg',trleg,'leglocation','northeast');
        xlabel('Block number'); ylabel('Peak press force (N)'); set(gca,'fontsize',fs); axis square;
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
        
    case 'press_dur' % analysis of press duration and press overlap (dissecting IPIs into press2release and release2press)
        if nargin > 1 % load single subj data
            subj = varargin{1};
            D = load( fullfile(pathToData, sprintf('sr2_%s.mat',subj)) ); %load data for this subject
            D.SN = ones( numel(D.TN), 1) * str2double(varargin{1}(2:3)); %add information about subject number
        else % load group data
            D = load( fullfile(pathToAnalyze, 'sr2_training_all_data.mat') );
        end
        
        % create summary table
        T = tapply(D,{'SN', 'day', 'BN', 'train'},...
            {D.rel2press_dur>0,'sum', 'name','gapsum'},...
            {D.rel2press_dur<5000,'sum', 'name','gaptot'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        T.gapper = T.gapsum./T.gaptot;
        
        % open figure
        if nargin>1; figure('Name',sprintf('Gap percentage - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Gap percentage - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        %T = normData(T, {'gapper'}, 'sub');
        
        subplot(2,2,1); title('');
        plt.line([T.day T.BN], mean(T.gapper,2)*100, 'plotfcn','nanmean', 'split',T.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg);
        xlabel('Block number'); ylabel('Percentage of transitions with gap'); set(gca,'fontsize',fs); axis square;
        
        subplot(2,2,2); title('');
        plt.line([T.day T.BN], T.gapper*100, 'plotfcn','nanmean', 'errorbars','shade', 'style',trnssty, 'leg',{'T1','T2','T3','T4'});
        xlabel('Block number'); ylabel('Percentage of transitions with gap'); set(gca,'fontsize',fs); axis square;
        
        D.rel2press_dur(D.rel2press_dur<=0) = NaN; % remove overlap transitions
        % create summary table
        T = tapply(D,{'SN', 'train', 'day', 'BN'},...
            {D.press2rel_dur,'nanmean', 'name','p2r'},...
            {D.rel2press_dur,'nanmean', 'name','r2p'},...
            'subset', D.isError==0 & D.dummy==0 & D.rtt==0);
        
        % open figure
        if nargin>1; figure('Name',sprintf('Press duration / overlap - subj %02d',str2double(varargin{1}(2:3)))); else; figure('Name',sprintf('Press duration / overlap - group (N=%d)',ns)); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        T1.p2r           = reshape(T.p2r, [], 1);
        T1.SN            = repmat(T.SN, numel(T.p2r(1,:)), 1);
        T1.BN            = repmat(T.BN, numel(T.p2r(1,:)), 1);
        T1.day           = repmat(T.day, numel(T.p2r(1,:)), 1);
        T1.train         = repmat(T.train, numel(T.p2r(1,:)), 1);
        T1.T             = reshape(repmat(1:numel(T.p2r(1,:)), numel(T.p2r(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T1 = normData(T1, {'p2r'}, 'sub');
        
        subplot(2,2,1); title('');
        plt.line([T1.day T1.BN], T1.normp2r, 'plotfcn','nanmean', 'split',T1.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg);
        xlabel('Block number'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        subplot(2,2,2); title('');
        plt.line([T1.day T1.BN], T1.normp2r, 'plotfcn','nanmean', 'split',T1.T, 'errorbars','shade', 'style',presssty, 'leg',{'P1','P2','P3','P4','P5'});
        xlabel('Block number'); ylabel('Press duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        T2.r2p           = reshape(T.r2p, [], 1);
        T2.SN            = repmat(T.SN, numel(T.r2p(1,:)), 1);
        T2.BN            = repmat(T.BN, numel(T.r2p(1,:)), 1);
        T2.day           = repmat(T.day, numel(T.r2p(1,:)), 1);
        T2.train         = repmat(T.train, numel(T.r2p(1,:)), 1);
        T2.T             = reshape(repmat(1:numel(T.r2p(1,:)), numel(T.r2p(:,1)), 1), [], 1);
        
        % normalize MT data to remove between-subject variability (i.e. plot within-subject standard error)
        T2 = normData(T2, {'r2p'}, 'sub');
        
        subplot(2,2,3); title('');
        plt.line([T2.day T2.BN], T2.normr2p, 'plotfcn','nanmean', 'split',T2.train, 'errorbars','shade', 'style',trsty_cbs, 'leg',trleg);
        xlabel('Block number'); ylabel('Gap duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        subplot(2,2,4); title('');
        plt.line([T2.day T2.BN], T2.normr2p, 'plotfcn','nanmean', 'split',T2.T, 'errorbars','shade', 'style',trnssty, 'leg',{'T1','T2','T3','T4'});
        xlabel('Block number'); ylabel('Gap duration (ms)'); set(gca,'fontsize',fs); axis square;
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end
end

function [theta_hat]=fitlog(x,y,theta_zero)
fcn=@(theta)bernoulli_loss(y,modlog(theta,x));
theta_hat=fminsearch(fcn,theta_zero);
end

function [y_hat]=modlog(theta,x)
a=theta(1,1); b=theta(2,1);
%c=theta(3,1);
c=.20;
y_hat=1./(1 + exp(-a*x + b))*(1-c) + c;
end

function [loss]=bernoulli_loss(y,y_hat)
n=numel(y);
loss=zeros(n,1);
for i=1:n; loss(i,1)=(y(i)*log(y_hat(i)) + (1-y(i))*log(1-y_hat(i))); end
loss=-(nansum(loss));
end