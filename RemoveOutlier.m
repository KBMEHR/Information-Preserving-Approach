function [L,R] = RemoveOutlier(L,R)

intLeng = R-L;
left = sort(L);
right = sort(R);
leng = sort(intLeng);
n=length(L);

NN1 = floor(n * 0.25 + 0.5);
NN2 = floor(n * 0.75 + 0.5);

% Compute Q(0.25), Q(0.75) and IQR for left endpoints
QL25 = (0.5 - n * 0.25 + NN1) * left(NN1) + (n * 0.25 + 0.5 - NN1) * left(NN1+1);
QL75 = (0.5 - n * 0.75 + NN2) * left(NN2) + (n * 0.75 + 0.5 - NN2) * left(NN2+1);
LIQR = QL75 - QL25;

% Compute Q(0.25), Q(0.75) and IQR for right-endpoints
QR25 = (0.5 - n * 0.25 + NN1) * right(NN1) + (n * 0.25 + 0.5 - NN1) * right(NN1+1);
QR75 = (0.5 - n * 0.75 + NN2) * right(NN2) + (n * 0.75 + 0.5 - NN2) * right(NN2+1);
RIQR = QR75 - QR25;

% Compute Q(0.25), Q(0.75) and IQR for interval length
QLeng25 = (0.5 - n * 0.25 + NN1) * leng(NN1) + (n * 0.25 + 0.5 - NN1) * leng(NN1+1);
QLeng75 = (0.5 - n * 0.75 + NN2) * leng(NN2) + (n * 0.75 + 0.5 - NN2) * leng(NN2+1);
lengIQR = QLeng75 - QLeng25;
bound=.25;

% Left & Right outlier processing
for i=n:-1:1
    if (LIQR>bound & (L(i)<QL25-1.5*LIQR | L(i)>QL75+1.5*LIQR))...
            |(RIQR>bound & (R(i)<QR25-1.5*RIQR | R(i)>QR75+1.5*RIQR))
        L(i) = [];
        R(i) = [];
        intLeng(i)=[];
    end
end

% numSurvivedData = length(L);
end

