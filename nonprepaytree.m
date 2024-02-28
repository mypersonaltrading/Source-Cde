function output=nonprepaytree(fv,rmgte,n,freq,shortTree)

r=(1+rmgte)^(1/freq)-1;           
T=n*freq;                          
dt=1/freq;                     
non_prepay_table=nan(T+1,T+1);      
mgtepmt=fv*r/(1-1/(1+r)^T);    

non_prepay_table(T+1,:)=mgtepmt;      
for i=T:-1:1
    for j=1:i
        if i==1
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt);
        else
            non_prepay_table(i,j)=0.5*(non_prepay_table(i+1,j)+non_prepay_table(i+1,j+1))*exp(-shortTree(i,j)*dt)+mgtepmt;
        end        
    end
end



p0=non_prepay_table(1,1);      
pU=non_prepay_table(2,2);     
pD=non_prepay_table(2,1);   
pUU=non_prepay_table(3,3);
pUD=non_prepay_table(3,2);
pDD=non_prepay_table(3,1);

yUU=shortTree(3,3);        
yUD=shortTree(3,2);
yDD=shortTree(3,1);
yU=shortTree(2,2) ;
yD=shortTree(2,1);

duration=-1/p0*(pU-pD) /(yU-yD);
duraUp=-1/pU*(pUU-pUD) / (yUU-yUD);
duraDown=-1/pD*(pUD-pDD)/(yUD-yDD);
convexity=duration^2-(duraUp-duraDown)/(yU-yD);
output=[p0,duration,convexity];
%output=non_prepay_table;