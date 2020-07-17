function mockBarPlots(mTarget,stdTarget,mNon,stdNon,n,YLIM)

%Call up figure and format
figure1 = figure;
axes1 = axes('Parent',figure1,'YTickLabel',{'0','0.5','1'},...
    'YTick',[0 .5 1.0],...
    'XTick',zeros(1,0),'XLim',[-1,2],'YLim',YLIM);

hold on

%Generate target starts
targets = (randn(n,1)*stdTarget)+mTarget;

%Generate non-target starts
nons = (randn(n,1)*stdNon)+mNon;

%Plot the mean and std bars
bar([0],[mTarget],'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
bar([1],[mNon],'FaceColor','none','EdgeColor',[0 0 0],'BarWidth',0.6);
errorbar(0,mTarget,stdTarget,'Marker','.','LineStyle','none','Color',[0 0 0])
errorbar(1,mNon,stdNon,'Marker','.','LineStyle','none','Color',[0 0 0])

%Plot the individual points on top
for i = 1:n
    plot(0,targets(i),'ok')
    plot(1,nons(i),'ok')
end