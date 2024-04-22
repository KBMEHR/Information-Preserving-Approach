clc; clear all; close all

%% Read Data
[A,words] = xlsread('OrderedWEBdata.xls');
words=words(1,1:2:end);
[row, col] = size(A);
L = zeros(row,col/2);
R = zeros(row,col/2);

%% Determine optimal MFs (non-parametric)
step = 0.01;
elements = 0:step:10; 
DiscrtMFs = zeros(numel(words),numel(elements),3);
for i = 1:length(words)
	L(:,i) = A(1:row, 2*i-1);   % Left end-points for interval data
    R(:,i) = A(1:row, 2*i);     % Right end-points for interval data
    DiscrtMFs(i,:,:) = DetermineFOU(L(:,i), R(:,i), elements);
end

%% Plot non-parametric FOUs
figure
for i=1:length(words)
    subplot(8,4,i);
    hold on
    plot([0 0 10 10 0],[0 1 1 0 0],'k')
    plot(elements,DiscrtMFs(i,:,1),'b')
    plot(elements,DiscrtMFs(i,:,3),'g')
    plot(elements,DiscrtMFs(i,:,2),'r')
    title(words(i),'fontsize',8)
    set(gca,'XTick',[],'YTick',[])
end

%% Approximate non-parametric FOUs by trapezoidal FOUs
words = convertCharsToStrings(words);
disp("Trapezoidal FOU for words:")
TrapezoidalFOU = zeros(length(words),9);
Cs = zeros(length(words),1);
for i = 1:length(words)
    MedianMFs = DiscrtMFs(i,:,3);
    Lmed = elements(find(MedianMFs>0,1,'first'));
    Rmed = elements(find(MedianMFs>0,1,'last'));
    medRange = Rmed-Lmed;
    UMF_LowerBound = [0, 0, 0, 0];  
    UMF_UpperBound = [Lmed,10,10,10];   
    LMF_LowerBound = [Lmed, 0, 0, 0];  
    LMF_UpperBound = [Rmed, medRange, medRange, medRange];   
    MFs = [elements', DiscrtMFs(i,:,1)', DiscrtMFs(i,:,2)'];
    TrapezoidalFOU(i,:) = ApproxTrapFOU(MFs, UMF_LowerBound, UMF_UpperBound, LMF_LowerBound, LMF_UpperBound);
    Cs(i)=centroidIT2(TrapezoidalFOU(i,:));
    disp(strcat(words(i),": ",num2str(TrapezoidalFOU(i,:))))
end
[Cs,index]=sort(Cs);    % Sort the centers of the centroids
TrapezoidalFOU = TrapezoidalFOU(index,:);     
words = words(index);

%% Plot final FOUs
figure
for i=1:length(words)
    subplot(8,4,i);
    fill(TrapezoidalFOU(i,[1:4 8:-1:5 1]),[0 1 1 0 0 TrapezoidalFOU(i,[9 9]) 0 0],[.9 .9 .9]) 
	title(strcat(words(i),''),'fontsize',8)
    set(gca,'XTick',[],'YTick',[])
    axis([0 10 0 1])
end
