analysisfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/AnalysisFigs/CTXSlowOscillationAnalysis';
SlowOscillationAll = GetMatResults(analysisfolder,'_CTXSlowOscillationAnalysis');
figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/AnalysisFigs/CTXSlowOscillationAnalysis';
numrecs = length(SlowOscillationAll);
SlowOscillationAll = CollapseStruct(SlowOscillationAll);
%%

DOWNpeakLFP = CollapseStruct(SlowOscillationAll.DOWNpeakLFP,2,'justcat',true);
DOWNpeakLFP.t = SlowOscillationAll.DOWNpeakLFP.t;

celldeltaPETH = CollapseStruct(SlowOscillationAll.celldeltaPETH,2,'justcat',true);
[~,celldeltaPETH.sortE] = sort(mean(celldeltaPETH.pE,1));
[~,celldeltaPETH.sortI] = sort(mean(celldeltaPETH.pI,1));
celldeltaPETH.numE = length(celldeltaPETH.sortE);
celldeltaPETH.numI = length(celldeltaPETH.sortI);

celldeltaPETH.ALL = [celldeltaPETH.pE celldeltaPETH.pI];
[~,celldeltaPETH.sortALL] = sort(mean(celldeltaPETH.ALL,1));
celldeltaPETH.numALL = length(celldeltaPETH.sortALL);

popdeltaPETH = CollapseStruct(SlowOscillationAll.popdeltaPETH,2,'justcat',true);

UPDOWNratehist = CollapseStruct(SlowOscillationAll.UPDOWNratehist,1,'justcat',true);

UPDOWNdurhist = CollapseStruct(SlowOscillationAll.UPDOWNdurhist,1,'justcat',true);

specgram = CollapseStruct(SlowOscillationAll.meanFFT,2,'justcat',true);
specgram.freqs=SlowOscillationAll.meanFFT.freqs;

LFPhist = CollapseStruct(SlowOscillationAll.LFPhist,1);

%% UP/DOWN Stats


%%
Ecolor = makeColorMap([0 0 0],[0 0.5 0],[1 1 1]);
Icolor = makeColorMap([0 0 0],[0.7 0 0],[1 1 1]);
ALLcolor = makeColorMap([1 1 1],[0.5 0.5 0.5],[0 0 0]);
xwin = [-0.7 0.7];

UPcolor = [0.7 0.5 0];
UPcolor = [0.7 0 0];
DOWNcolor = [0 0 0.8];

figure
    subplot(6,2,1)
        plot(DOWNpeakLFP.t,mean(DOWNpeakLFP.mean,2),'k','linewidth',2)
        hold on
        errorshade(DOWNpeakLFP.t,mean(DOWNpeakLFP.mean,2),std(DOWNpeakLFP.mean,[],2),std(DOWNpeakLFP.mean,[],2),'k','scalar')        
        axis tight
        plot([0 0],get(gca,'ylim'),'b--','linewidth',0.1)
        box off
        xlim(xwin)
        set(gca,'xtick',[]);set(gca,'ytick',[])
        ylabel('LFP')
        
% subplot(3,2,3)
%     imagesc(celldeltaPETH.t(:,1),[0 1],log2(celldeltaPETH.pE(:,celldeltaPETH.sortE))')
%     hold on
%     axis xy
%     plot([0 0],get(gca,'ylim'),'w')
%     colormap(gca,Ecolor)
%     xlim(xwin)
%     %colorbar
%     set(gca,'xtick',[]);set(gca,'ytick',[])
%     ylabel(['pE Cells (',num2str(celldeltaPETH.numE),')'])
% subplot(6,2,3)
%     imagesc(celldeltaPETH.t(:,1),[0 1],log2(celldeltaPETH.pI(:,celldeltaPETH.sortI))')
%     hold on
%     axis xy
%     plot([0 0],get(gca,'ylim'),'w')
%     xlim(xwin)
%     colormap(gca,Icolor)
%     %colorbar
%     set(gca,'xtick',[]);set(gca,'ytick',[])
%     ylabel(['pI Cells (',num2str(celldeltaPETH.numI),')'])


subplot(6,2,[3 5])
    imagesc(celldeltaPETH.t(:,1),[0 1],log10(celldeltaPETH.ALL(:,celldeltaPETH.sortALL))')
    hold on
    axis xy
    plot([0 0],get(gca,'ylim'),'b--','linewidth',0.1)
    xlim(xwin)
    colormap(gca,ALLcolor)
    colorbar
    caxis([-1.5 1])
    set(gca,'xtick',[]);set(gca,'ytick',[])
    ylabel(['pI Cells (',num2str(celldeltaPETH.numALL),')'])


% subplot(4,2,7)
%     plot(popdeltaPETH.t(:,1),4.*nanmean(popdeltaPETH.pE,2),'color',Ecolor(end/2,:),'linewidth',2)
% 	hold on
%     plot(popdeltaPETH.t(:,1),nanmean(popdeltaPETH.pI,2),'color',Icolor(end/2,:),'linewidth',2)
%     xlim(xwin)
%     ylabel('Mean Pop. Rate (Hz)')
%     xlabel('t (s - aligned to delta peak)')
%     box off
%     axis tight
    
    
subplot(4,2,5)
    plot(popdeltaPETH.t(:,1),mean(popdeltaPETH.ALL,2),'color',[0.5 0.5 0.5],'linewidth',2)
	hold on
    errorshade(popdeltaPETH.t(:,1),mean(popdeltaPETH.ALL,2),std(popdeltaPETH.ALL,[],2),std(popdeltaPETH.ALL,[],2),'k','scalar')
    plot([0 0],get(gca,'ylim'),'b--','linewidth',0.1)
    xlim(xwin)
    ylabel('Mean Pop. Rate (Hz)')
    xlabel('t (s - aligned to delta peak)')
    box off
    axis tight
    
    
% subplot(3,2,2)
%     plot(UPDOWNratehist.interpbins,mean(UPDOWNratehist.interpNREM,1),'k','linewidth',2)
%     hold on
%     plot(UPDOWNratehist.interpbins,mean(UPDOWNratehist.interpDOWN,1),'color',DOWNcolor,'linewidth',1)
%     %hold on
%     plot(UPDOWNratehist.interpbins,mean(UPDOWNratehist.interpUP,1),'color',UPcolor,'linewidth',1)
%     box off
%     axis tight
    
    exhist = 8;
subplot(5,2,2)
%plot(UPDOWNratehist.interpbins(1,:),mean(UPDOWNratehist.interpNREM(numcells>60,:),1))
% bar(UPDOWNratehist.interpbins(1,:),(UPDOWNratehist.interpNREM(exhist,:)),'facecolor',[0.5 0.5 0.5])
%     hold on
%     plot(UPDOWNratehist.interpbins,(UPDOWNratehist.interpDOWN(exhist,:)),'color',DOWNcolor,'linewidth',2)
%     %hold on
%     plot(UPDOWNratehist.interpbins,(UPDOWNratehist.interpUP(exhist,:)),'color',UPcolor,'linewidth',2)
%   
bar(SlowOscillationAll.UPDOWNratehist(exhist).bins,(SlowOscillationAll.UPDOWNratehist(exhist).NREM),'facecolor',[0.5 0.5 0.5])
    hold on
    plot(SlowOscillationAll.UPDOWNratehist(exhist).bins,(SlowOscillationAll.UPDOWNratehist(exhist).DOWN),'--','color',DOWNcolor,'linewidth',3)
    %hold on
    plot(SlowOscillationAll.UPDOWNratehist(exhist).bins,(SlowOscillationAll.UPDOWNratehist(exhist).UP),'--','color',UPcolor,'linewidth',3)
    
    box off
    axis tight
    ylim([0 0.18])
   exhist = 1;
   
   
% subplot(3,2,4)
%     plot(UPDOWNdurhist.DOWNbins,UPDOWNdurhist.DOWNhist(exhist,:),'color',DOWNcolor,'linewidth',2)
%     hold on
%     plot(UPDOWNdurhist.UPbins,UPDOWNdurhist.UPhist(exhist,:),'color',UPcolor,'linewidth',2)
%     box off
%     axis tight
% 
%     box off
%     ylim([0 0.32])
%     %axis tight
    
subplot(3,2,4)
    plot(UPDOWNdurhist.logbins,mean(UPDOWNdurhist.DOWNlog,1),'color',DOWNcolor,'linewidth',2)
    hold on
    plot(UPDOWNdurhist.logbins,mean(UPDOWNdurhist.UPlog,1),'color',UPcolor,'linewidth',2)
    errorshade(UPDOWNdurhist.logbins(1,:),mean(UPDOWNdurhist.UPlog,1),std(UPDOWNdurhist.UPlog,[],1),std(UPDOWNdurhist.UPlog,[],1),UPcolor,'scalar')
    errorshade(UPDOWNdurhist.logbins(1,:),mean(UPDOWNdurhist.DOWNlog,1),std(UPDOWNdurhist.DOWNlog,[],1),std(UPDOWNdurhist.DOWNlog,[],1),DOWNcolor,'scalar')

    box on
    ylim([0 0.25])
    axis tight
    ylim([0 0.25])
    xlim([-1.5 1.5])
    LogScale('x',10)
    
    
 subplot(3,3,9)

    hold on
    plot(log2(specgram.freqs),mean(specgram.REM,2),'r','linewidth',2)
    plot(log2(specgram.freqs),mean(specgram.WAKE,2),'k','linewidth',2)
    plot(log2(specgram.freqs),mean(specgram.NREM,2),'b','linewidth',2)
    
    errorshade(log2(specgram.freqs),mean(specgram.REM,2),std(specgram.REM,[],2),std(specgram.REM,[],2),'r','scalar')
    errorshade(log2(specgram.freqs),mean(specgram.WAKE,2),std(specgram.WAKE,[],2),std(specgram.WAKE,[],2),'k','scalar')
    plot(log2(specgram.freqs),mean(specgram.NREM,2),'b','linewidth',2)
    errorshade(log2(specgram.freqs),mean(specgram.NREM,2),std(specgram.NREM,[],2),std(specgram.NREM,[],2),'b','scalar')
    axis tight
    box off
    LogScale('x',2)
    xlim([log2(0.25) log2(100)])
    
    
subplot(3,3,8)
    plot(LFPhist.bins(1,:),mean(LFPhist.REM,1),'r','linewidth',2)
    hold on
    plot(LFPhist.bins(1,:),mean(LFPhist.WAKE,1),'k','linewidth',2)
    plot(LFPhist.bins(1,:),mean(LFPhist.NREM,1),'b','linewidth',2)
    
    errorshade(LFPhist.bins(1,:),mean(LFPhist.REM,1),std(LFPhist.REM,[],1),std(LFPhist.REM,[],1),'r','scalar')
    errorshade(LFPhist.bins(1,:),mean(LFPhist.WAKE,1),std(LFPhist.WAKE,[],1),std(LFPhist.WAKE,[],1),'k','scalar')
    plot(LFPhist.bins(1,:),mean(LFPhist.NREM,1),'b','linewidth',2)
    errorshade(LFPhist.bins(1,:),mean(LFPhist.NREM,1),std(LFPhist.NREM,[],1),std(LFPhist.NREM,[],1),'b','scalar')
    xlim([-0.25e4 0.25e4])
    ylim([0 0.07])
    box off

    
NiceSave('CTXSlowOscillation',analysisfolder,[])


%%
for rr = 1:numrecs
    numcells(rr) = length(SlowOscillationAll.celldeltaPETH(rr).pE(1,:))+length(SlowOscillationAll.celldeltaPETH(rr).pI(1,:));
end
[~,sortnumcells] = sort(numcells);
%for ee =1:numrecs
exhist = 8;%exhist = 15;
figure
subplot(2,2,1)
plot(UPDOWNratehist.interpbins(exhist,:),UPDOWNratehist.interpNREM(exhist,:),'k','linewidth',2)

subplot(2,2,2)
imagesc(UPDOWNratehist.interpNREM(sortnumcells,:))

subplot(2,2,3)
%plot(UPDOWNratehist.interpbins(1,:),mean(UPDOWNratehist.interpNREM(numcells>60,:),1))
    plot(UPDOWNratehist.interpbins,(UPDOWNratehist.interpNREM(exhist,:)),'k','linewidth',2)
    hold on
    plot(UPDOWNratehist.interpbins,(UPDOWNratehist.interpDOWN(exhist,:)),'color',DOWNcolor,'linewidth',1)
    %hold on
    plot(UPDOWNratehist.interpbins,(UPDOWNratehist.interpUP(exhist,:)),'color',UPcolor,'linewidth',1)
    box off
    axis tight

subplot(2,2,4)
for nn = 10:10:100
hold on
plot(UPDOWNratehist.interpbins(1,:),mean(UPDOWNratehist.interpNREM(numcells>nn,:),1))
end
%end