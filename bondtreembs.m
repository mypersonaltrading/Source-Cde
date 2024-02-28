function modpandvol=bondtreembs(shortTree,N,startParam,dt)
%function [modPrice,modVol]=bondtree(shortTree,T,startParam)
BT=ones(N+1,N+1); 

sT=shortTree;
%for i=1:N
    for j=1:N
        if j==1
            sT(N,j)=startParam(1);
        else
            sT(N,j)=sT(N,j-1)*exp(2*startParam(2)*sqrt(dt));
        end
    end
%end

for i=N:-1:1
    for j=1:i
        BT(i,j)=(0.5*BT(i+1,j)+0.5*BT(i+1,j+1))*exp(-sT(i,j)*dt);
    end
end
%output for the bondtree function
modPrice=BT(1,1);
upVol=-(1/(N-dt))*log(BT(2,2));
downVol=-(1/(N-dt))*log(BT(2,1));
modVol=0.5/sqrt(dt)*log(upVol/downVol);
modpandvol=[modPrice,modVol];






