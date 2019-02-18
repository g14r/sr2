function [D]=sr2_training_trial(MOV,D,fig,fig_name,varargin)
%% function [D]=sr2_training_trial(MOV,D,fig,fig_name,varargin)
% Trial routine for SequenceRepetition experiment (sr2_training), training variant
% called by the Subject routine (sr2_training_subj.m)

%% extract data
if (isempty(MOV))
    return;
end;
%state=MOV(:,1);
%absTime=MOV(:,2);
%time=MOV(:,3);
time=2:2:(size(MOV(:,1),1))*2;
%force=smooth_kernel(MOV(:,4:end),4);
force=MOV(:,4:end);
RH=[6 7 8 9 10]; %right hand column indices
LH=[1 2 3 4 5]; %left hand column indices

%%
fNames=fieldnames(D);
resp=[]; %which presses?

%%
for i=1:length(fNames)
    if length(fNames{i})>8
        if strcmp(fNames{i}(1:8),'response')
            eval(['resp=[resp,D.',fNames{i},'];']);
        end
    end
end
D.numPress=sum(resp>0); %number of presses
if D.numPress==5
    pressTime=zeros(1,D.numPress);
    relTime=zeros(1,D.numPress);
    for press=1:D.numPress
        pressNum=['D.pressTime',num2str(press)];
        pressTime(press)=eval(pressNum); %time of press
        relTime(press)=time(find(force(:,RH(eval(['D.response',num2str(press)])))>=.95, 1, 'last')); %time of release
        eval(['D.relTime',num2str(press),'=',num2str(relTime(press)),';']);
    end
    D.pressTimes=pressTime;
    D.relTimes=relTime;
    D.IPI=diff(pressTime);
    D.press2rel_dur=relTime-pressTime; %duration from press to press release
    D.rel2press_dur=pressTime(2:end)-relTime(1:end-1); %duration from press release to subsequent press (negative values mean that subsequent press happened before previous release; positive values vice versa)
else
    pressTime=zeros(1,D.numPress);
    relTime=zeros(1,D.numPress);
    for press=1:D.numPress
        pressNum=['D.pressTime',num2str(press)];
        pressTime(press)=eval(pressNum); %time of press
        relTime(press)=time(find(force(:,RH(eval(['D.response',num2str(press)])))>=.95, 1, 'last')); %time of release
        eval(['D.relTime',num2str(press),'=',num2str(relTime(press)),';']);
    end
end

%% Display trial
if (fig>0)
    figure('Name',fig_name)
    %subplot(2,1,1)
    if D.hand==1
        plot(time,force(:,LH),'LineWidth',2);
        title('Force traces for LEFT hand presses','FontSize',20);
    elseif D.hand==2
        plot(time,force(:,RH),'LineWidth',2);
        title('Force traces for RIGHT hand presses','FontSize',20);
    end
    xlabel('Time (ms)'); ylabel('Force (N)'); set(gca,'FontSize',20); xlim([2500 4500]); ylim([-0.5 4.5]); axis square;
    hold on;
    drawline(pressTime,'dir','vert', 'linestyle','-', 'color','b');
    drawline(relTime,'dir','vert', 'linestyle','-', 'color','r');
    drawline(3100,'dir','vert', 'linestyle','--');
    drawline(3300,'dir','vert', 'linestyle','--');
    drawline(1,'dir','horz', 'linestyle','--');
    legend({'Thumb','Index','Middle','Ring','Little'},'FontSize',20)
    hold off;
    %     subplot(2,1,2)
    %     plot(time,state,'LineWidth',2);
    %     ylim([1 7]);
    %     xlabel('Time (ms)'); ylabel('State'); set(gca,'FontSize',20);
    %     hold on;
    %     drawline(pressTime,'dir','vert'); legend state
    %     KbWait;
end