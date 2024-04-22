function MFs = DetermineFOU(L, R, points)
% L: left end-points of the intervals from survey
% R: right end-points of the intervals from survey
% points: elements in the discritized domain
% MFs: MFs of Optimal (non-parametric) LMF and UMF along with Median MF

%% Remove invalid data 
index=find(isnan(L)+isnan(R)); % remove incomplete data
L(index)=[];
R(index)=[];
for i=length(L):-1:1
    if L(i)<0 | L(i)>10 | R(i)<0 | R(i)>10 |  R(i)<=L(i) | R(i)-L(i)>=10
        L(i) = [];
        R(i) = [];
    end
end

%% Outlier processing based on left & right endpoints
intLeng = R-L;
left = sort(L);
right = sort(R);
n=length(L);

NN1 = floor(n * 0.25 + 0.5);
NN2 = floor(n * 0.75 + 0.5);
QL25 = (0.5 - n * 0.25 + NN1) * left(NN1) + (n * 0.25 + 0.5 - NN1) * left(NN1+1);
QL75 = (0.5 - n * 0.75 + NN2) * left(NN2) + (n * 0.75 + 0.5 - NN2) * left(NN2+1);
LIQR = QL75 - QL25;
QR25 = (0.5 - n * 0.25 + NN1) * right(NN1) + (n * 0.25 + 0.5 - NN1) * right(NN1+1);
QR75 = (0.5 - n * 0.75 + NN2) * right(NN2) + (n * 0.75 + 0.5 - NN2) * right(NN2+1);
RIQR = QR75 - QR25;

bound=.25;
for i=n:-1:1
    if (LIQR>bound & (L(i)<QL25-1.5*LIQR | L(i)>QL75+1.5*LIQR))...
            |(RIQR>bound & (R(i)<QR25-1.5*RIQR | R(i)>QR75+1.5*RIQR))
        L(i) = [];
        R(i) = [];
        intLeng(i)=[];
    end
end
n = length(L);

%% Map data intervals to triangular fuzzy sets 
T1L = 0.5*(L+R) - sqrt(2)*(R-L)/2;
T1R = 0.5*(L+R) + sqrt(2)*(R-L)/2;
T1M = (T1L+T1R)./2;

% figure
% hold on
% axis([0 10 0 1])
T1MFs = zeros(n,numel(points));
for i =1:n
    T1MFs(i,:) = trimf(points, [T1L(i),T1M(i),T1R(i)]);   
%     plot([T1L(i),T1M(i),T1R(i)],[0 1 0],'color',rand(1,3),'linewidth',1);
end

%% Determine optimal FOUs 
MFs = zeros(numel(points),3);
for i = 1:numel(points)
    pointMFs = sort(T1MFs(:,i));
    pointMedian = median(pointMFs);
    numLowerMFs = n/2;
    if mod(n,2) == 0    % Number of MFs is even
        pointLowerMFs = pointMFs(1:numLowerMFs);  
        pointUpperMFs = pointMFs(numLowerMFs+1:n);  
    else
        numLowerMFs = floor(numLowerMFs);
        pointLowerMFs = pointMFs(1:numLowerMFs);
        pointUpperMFs = pointMFs(numLowerMFs+2:n);  
    end    
    numUpperMFs = numLowerMFs;
    pointLowerMFs = flipud(pointLowerMFs);  
        
    MFsInds = 1:numUpperMFs;
    Coverage = (MFsInds./numUpperMFs)';
   
    UpperMFsDisp = pointUpperMFs-pointMedian;
    if any(UpperMFsDisp)
      UpperSpecificity = 1 - UpperMFsDisp./max(UpperMFsDisp);
    else
       UpperSpecificity(1:length(UpperMFsDisp),1) = 1;
    end
    UpperPerformance = Coverage.*UpperSpecificity;
    [~,indUMF] = max(UpperPerformance);
    pointUMF = pointUpperMFs(indUMF);

    LowerMFsDisp = pointMedian-pointLowerMFs;
    if any(LowerMFsDisp)
       LowerSpecificity = 1 - LowerMFsDisp./max(LowerMFsDisp);
    else
       LowerSpecificity(1:length(LowerMFsDisp),1) = 1;
    end
    LowerPerformance = Coverage.*LowerSpecificity ;
    [~,indLMF] = max(LowerPerformance);
    pointLMF = pointLowerMFs(indLMF);

    MFs(i,:) = [pointLMF, pointUMF, pointMedian];
end

end

