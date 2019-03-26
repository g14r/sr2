function [xmodel,ymodel]=sr2_srttPilot_analyze
% function [xmodel,ymodel]=sr2_srttPilot_analyze
%
% example call: [xmodel,ymodel]=sr2_srttPilot_analyze;
%%
pathToData='/Users/gariani/Documents/data/SequenceRepetition/sr2';
%svec=[1,2,3,4,5,6,7,8];
svec=[1,2,3,4,5,6,7,8];
fs=16;
D=struct();
for s=svec
    fn=sprintf('sr2_srt_s%02d.dat',s);
    fname=fullfile(pathToData,fn);
    S=dload(fname);
    S.SN=ones(numel(S.TN),1)*s;
    D=addstruct(D,S);
end
%%
D1=tapply(D,{'SN','prepTime'},...
    {D.timingError,'mean','name','mTE'},...
    'subset',D.dummy==0 & D.prepTime>=100);
figure;
subplot(1,2,1); title(sprintf('Proportion of wrong timing trials (regardless of finger selection), N=%d',numel(svec)))
lineplot(D1.prepTime,(D1.mTE)*100,'errorbars','shade','style_thickline');
drawline(0,'dir','horz','linestyle','--'); drawline(25,'dir','horz','linestyle','--'); drawline(50,'dir','horz','linestyle','--'); ylim([-10 110]);drawline(75,'dir','horz','linestyle','--'); drawline(100,'dir','horz','linestyle','--');
drawline(100,'dir','vert','linestyle','--');drawline(200,'dir','vert','linestyle','--');drawline(300,'dir','vert','linestyle','--');drawline(400,'dir','vert','linestyle','--');drawline(500,'dir','vert','linestyle','--'); xlim([100 2400]); axis square; 
xlabel('Preparation time (ms)'); ylabel('Timing error (%)'); set(gca,'fontsize',fs); xticks([100:100:500,800,1600,2400]); yticks(0:10:100); legend('data','se');
%
subplot(1,2,2); title('Individual subjects')
lineplot(D1.prepTime,(D1.mTE)*100,'split',D1.SN,'errorbars','plusminus','errorcap',0,'style_thickline','leg','auto');
drawline(0,'dir','horz','linestyle','--'); drawline(25,'dir','horz','linestyle','--'); drawline(50,'dir','horz','linestyle','--'); ylim([-10 110]);drawline(75,'dir','horz','linestyle','--'); drawline(100,'dir','horz','linestyle','--');
drawline(100,'dir','vert','linestyle','--');drawline(200,'dir','vert','linestyle','--');drawline(300,'dir','vert','linestyle','--');drawline(400,'dir','vert','linestyle','--');drawline(500,'dir','vert','linestyle','--'); axis square; xlim([100 1600]);
xlabel('Preparation time (ms)'); ylabel('Timing error (%)'); set(gca,'fontsize',fs); xticks([100:100:500,800,1600,2400]); yticks(0:10:100);
%%
D2=tapply(D,{'SN','prepTime'},...
    {D.isError,'mean','name','mTE'},...
    'subset',D.dummy==0 & D.timingError==0 & D.prepTime>=100);
figure;
subplot(1,2,1); title(sprintf('Proportion of correct finger selection (within correct timing only), N=%d',numel(svec)))
[x,y]=lineplot(D2.prepTime,(1-D2.mTE),'errorbars','shade','style_thickline');
%% fit logistic function
theta_zero=[.02,250,0.25]';
[theta_hat]=fitlog(x,y,theta_zero);
[y_hat]=modlog(theta_hat,x);
hold on; [xmodel,ymodel]=lineplot(x,y_hat,'linecolor','r','linewidth',2,'errorbars','');
%
drawline(0,'dir','horz','linestyle','--'); drawline(.25,'dir','horz','linestyle','--'); drawline(.50,'dir','horz','linestyle','--'); ylim([-.10 1.10]);drawline(.75,'dir','horz','linestyle','--'); drawline(1.00,'dir','horz','linestyle','--');
drawline(.90,'dir','horz','linestyle','-','color','r'); drawline(.95,'dir','horz','linestyle','-','color','r');
drawline(100,'dir','vert','linestyle','--');drawline(200,'dir','vert','linestyle','--');drawline(300,'dir','vert','linestyle','--');drawline(400,'dir','vert','linestyle','--');drawline(500,'dir','vert','linestyle','--'); axis square; %xlim([100 2400]);  
xlabel('Preparation time (ms)'); ylabel('Selection accuracy (proportion)'); set(gca,'fontsize',fs); xticks([100:100:500,800,1600,2400]); yticks(0:.10:1.00); legend('data','','model');
%
subplot(1,2,2); title('Individual subjects')
lineplot(D2.prepTime,(1-D2.mTE)*100,'split',D2.SN,'errorbars','plusminus','errorcap',0,'style_thickline','leg','auto');
drawline(0,'dir','horz','linestyle','--'); drawline(25,'dir','horz','linestyle','--'); drawline(50,'dir','horz','linestyle','--'); ylim([-10 110]);drawline(75,'dir','horz','linestyle','--'); drawline(100,'dir','horz','linestyle','--');
drawline(100,'dir','vert','linestyle','--');drawline(200,'dir','vert','linestyle','--');drawline(300,'dir','vert','linestyle','--');drawline(400,'dir','vert','linestyle','--');drawline(500,'dir','vert','linestyle','--'); axis square; xlim([0 600]); 
xlabel('Preparation time (ms)'); ylabel('Selection accuracy (%)'); set(gca,'fontsize',fs); xticks([100:100:500,800,1600,2400]); yticks(0:10:100);
%%
% nurt=numel(unique(D.RT));
% figure; title('RT histogram (response window)');
% histplot(D.RT,'numcat',nurt,'percent',0); drawline(-100,'dir','vert','linestyle','--'); drawline(100,'dir','vert','linestyle','--'); xlim([-500 1000]); xlabel('RT'); ylabel('N trials'); set(gca,'fontsize',fs);
% figure;
% histplot(D.RT,'numcat',nurt,'percent',0,'splitrow',D.SN);
% figure;
% histplot(D.RT,'numcat',nurt,'percent',0,'splitrow',D.BN);
% figure;
% histplot(D.RT,'numcat',nurt,'percent',0,'splitrow',D.prepTime);
% figure; title('Relationship RT/preparation time');
% scatterplot(D.prepTime,D.RT,'regression','linear','markertype','o','markerfill','auto','markersize',9,'printcorr','subset',D.RT<=1000);
% xlabel('Preparation time'); ylabel('RT'); set(gca,'fontsize',fs); ylim([-500 1000]);
%%
end

function [theta_hat]=fitlog(x,y,theta_zero)
fcn=@(theta)bernoulli_loss(y,modlog(theta,x));
theta_hat=fminsearch(fcn,theta_zero);
end

function [y_hat]=modlog(theta,x)
a=theta(1,1); b=theta(2,1); c=theta(3,1);
y_hat=1./(1 + exp(-a*x + b))*(1-c) + c;
end

function [loss]=bernoulli_loss(y,y_hat)
n=numel(y);
loss=zeros(n,1);
for i=1:n
    loss(i,1)=(y(i)*log(y_hat(i)) + (1-y(i))*log(1-y_hat(i)));
end
loss=-(nansum(loss));
end
