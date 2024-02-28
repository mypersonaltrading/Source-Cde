clear;clc
%End date, in years
T=30;

%Number of steps in the binomial tree
N=360;

%set time step size
dt=T/N;

%import data for the treasury yield curve

yieldData=dataset('File','yieldAndVolatility2000.csv','Delimiter',',');

%par vector contains, respectively, beta0, beta1, beta2, and tau

startPar=[1,1,1,1];

% Define the Nelson-Siegel function
nsFunction = @(m, par) par(1) + (par(2) + par(3)) * ((1 - exp(-m/par(4))) ./ (m/par(4))) - par(3) * exp(-m/par(4));

%Alias the error function with observed data and parameters.

errorNelsonSiegelYield=@(par) sum((yieldData.y - nsFunction(yieldData.m, par)).^2);
errorNelsonSiegelVol=@(par) sum((yieldData.vol - nsFunction(yieldData.m, par)).^2);


%Run optimization
[parMinYield,errorMin] = fminsearch(errorNelsonSiegelYield, startPar);
[parMinVol,errorMin]   = fminsearch(errorNelsonSiegelVol, startPar);

%Creat a new set of maturities for monthly time steps over 30 years
monthlyMaturities=(dt:dt:T)';
modelleddata=dataset(monthlyMaturities);


%Alias new function fittedYield that takes as input some maturity monthlyMaturities and 
%returns fitted yields at the estimated parameters
fittedFunctionYield= @(monthlyMaturities) nsFunction(monthlyMaturities, parMinYield);
fittedFunctionVol  = @(monthlyMaturities) nsFunction(monthlyMaturities, parMinVol);

modelleddata.nelsonSiegelYield= fittedFunctionYield(monthlyMaturities);
modelleddata.nelsonSiegelVol= fittedFunctionVol(monthlyMaturities);




%backing out the price from the yield
modelleddata.nelsonSiegelPrice=exp(-modelleddata.nelsonSiegelYield.*modelleddata.monthlyMaturities);




%Create empty tree. the (i,j) node represents time i and node j of the
%short rate tree

shortTree=NaN(N,N);

%Starting Guesses for fminsearch
startR=.01;
startParam=[.01,.01];

i=1;
thisPrice=@(startR) exp(-startR*dt)*(0.5*1+0.5*1);
    thisMaturity=i-1;

%No volatility matching for first node

thisError=@(r) (thisPrice(r)-modelleddata.nelsonSiegelPrice(i))^2;

thisR=fminsearch(thisError,startR);

shortTree(i,1)=thisR;


%%Moving on in the tree
max_search=500*length(startParam);
options=optimset('MaxFunEvals',max_search,'MaxIter',max_search);

for i = 2:N
        i/N;
        thisPrice=modelleddata.nelsonSiegelPrice(i);
        thisMaturity=i;
        thisVolatility=modelleddata.nelsonSiegelVol(i);

        thisPriceAndVol=[thisPrice,thisVolatility];
        

        thisError=@(startParam) sum((thisPriceAndVol-bondtreembs(shortTree,i,startParam,dt)).^2);
        thisParam=fminsearch(thisError,startParam,options);
        startParam=thisParam;
        
        %Build next step in tree with the right parameters
        thisBottomR=thisParam(1);
        thisSigma=thisParam(2);

        treeStep=exp(2*thisSigma*sqrt(dt));

        %populate next step in shortTree
        shortTree(i,1)=thisBottomR; 
        for j=2:(i)
            shortTree(i,j)=shortTree(i,j-1)*treeStep;
        end

end



rmgte=0.07;
mrmgte=(1+rmgte)^(1/12)-1;  %monthly rate
fv=100000000;
freq=12;
mgtepmt=fv*(mrmgte*(1+mrmgte)^(360))/((1+mrmgte)^(360)-1);


%mortgage payment schedule
mgtepmtable=table(monthlyMaturities);
mgtepmtable.ipay(1)=fv*mrmgte;
mgtepmtable.ppalpay(1)=mgtepmt-mgtepmtable.ipay(1);
mgtepmtable.ppalout(1)=fv-mgtepmtable.ppalpay(1);

for t=2:N
    mgtepmtable.ipay(t)=mgtepmtable.ppalout(t-1)*mrmgte;
    mgtepmtable.ppalpay(t)=mgtepmt-mgtepmtable.ipay(t);
    mgtepmtable.ppalout(t)=mgtepmtable.ppalout(t-1)-mgtepmtable.ppalpay(t);
end



comp=mgtepmtable.ppalout(:);   %useful for prepayable tree buidling


%create the multiplication factors
factor=0.4:0.01:2.00;
%testing for factor 1
% factor=1;
% 
% test_npm=nonprepaytree(fv,rmgte,T,freq,shortTree);
% test_pm=prepaytree(fv,rmgte,T,freq,shortTree,comp);


for i =1:length(factor)
    prepayablemgte(i,:)=prepaytree(fv,rmgte,T,freq,shortTree.*factor(i),comp);          
    nonprepayablemgte(i,:)=nonprepaytree(fv,rmgte,T,freq,shortTree.*factor(i));
end



       
%% plotting 
subplot(3,2,1);
plot(yieldData.m, yieldData.y,'*',modelleddata.monthlyMaturities,modelleddata.nelsonSiegelYield,'LineWidth',2)
title('Fitted yield curve')
grid on
subplot(3,2,2);
grid on
plot(yieldData.m,yieldData.vol,'*',modelleddata.monthlyMaturities,modelleddata.nelsonSiegelVol,'LineWidth',2) 
title('Fitted volatility curve');
grid on

subplot(3,2,3)
plot(factor*shortTree(1,1),prepayablemgte(:,1),factor*shortTree(1,1),nonprepayablemgte(:,1), "LineWidth",2);
legend("Mtge value with prepayment option","Mtge value without prepayment option");
title("Value");
grid on

subplot(3,2,4)
plot(factor*shortTree(1,1),prepayablemgte(:,2),factor*shortTree(1,1),nonprepayablemgte(:,2), "LineWidth",2);
legend("Mtge duration with prepayment option","Mtge duration without prepayment option")
title("Duration");
grid on

subplot(3,2,5)
plot(factor*shortTree(1,1),prepayablemgte(:,3),factor*shortTree(1,1),nonprepayablemgte(:,3),'LineWidth',2);
legend("Mtge convexity with prepayment option","Mtge convexity without prepayment option")
title("Convexity");
grid on



