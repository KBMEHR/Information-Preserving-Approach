function ApproximatedFOU = ApproxTrapFOU(MFs, U_LB, U_UB, L_LB, L_UB)
% MFs: non-parametric FOU
% U_LB and U_UB: lower and upper bound for UMF in particle swarm optimization 
% L_LB and L_UB: lower and upper bound for LMF in particle swarm optimization 

clear -global UMF DiscrtUMF DiscrtLMF  
global elements DiscrtLMF DiscrtUMF UMF
elements = MFs(:,1);    % elements in the discritized domain
DiscrtLMF = MFs(:,2);   % Optimal (non-parametric) LMF 
DiscrtUMF = MFs(:,3);   % Optimal (non-parametric) UMF 

rng(5)
numParams = 4; 
options = optimoptions('particleswarm', 'FunctionTolerance',0.001, 'Display','off'); 
% options = optimoptions('particleswarm','MaxStallIterations',30,'FunctionTolerance',0.001, 'Display','off');   % slightly more accurate approximations
[OptParamUMF,ApprxErrUMF] = particleswarm(@FitnessUMF,numParams,U_LB,U_UB, options);
UMF = [OptParamUMF(1) sum(OptParamUMF(1:2)) sum(OptParamUMF(1:3)) sum(OptParamUMF(1:4))];
[OptParamLMF,ApprxErrLMF] = particleswarm(@FitnessLMF,numParams,L_LB,L_UB, options);
LMF = [OptParamLMF(1) sum(OptParamLMF(1:2)) sum(OptParamLMF(1:3)) sum(OptParamLMF(1:4)) max(DiscrtLMF)];
ApproximatedFOU = [UMF LMF];
end

function ApprxError = FitnessUMF(Param) 
  global elements DiscrtUMF 
  A = Param(1);
  B = A+Param(2);
  C = B+Param(3);
  D = C+Param(4);
  Penalty = 0;
  if D>10
    Penalty = max(1,(D-10)*10);
  end
  TrapMF = trapmf(elements,[A B C D]);
  ApprxError = mae(DiscrtUMF,TrapMF)+Penalty;       % Mean Absolute Error used as the criterion of approximation accuracy 
end

function ApprxError = FitnessLMF(Param) 
  global elements DiscrtLMF UMF
  A = Param(1);
  B = A+Param(2);
  C = B+Param(3);
  D = C+Param(4);
  maxMG = max(DiscrtLMF);
  Penalty = 0;
  if A < UMF(1)
    Penalty = Penalty+max(1,(UMF(1)-A)*10);
  end
  minValidB = UMF(1)+(UMF(2)-UMF(1))*maxMG;
  if B < minValidB
    Penalty = Penalty+max(1,(minValidB-B)*10);
  end
  maxValidC = UMF(4)-(UMF(4)-UMF(3))*maxMG;
  if C > maxValidC
    Penalty = Penalty+max(1,(C-maxValidC)*10);
  end
  if D > UMF(4)
      Penalty = Penalty+max(1,(D-UMF(4))*10);
  end
  TrapMF = trapmf(elements,[A B C D]);
  ApprxError = mae(DiscrtLMF,TrapMF*maxMG)+Penalty;    % Mean Absolute Error used as the criterion of approximation accuracy 
end

