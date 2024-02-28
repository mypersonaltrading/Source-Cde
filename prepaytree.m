function output=prepaytree(fv,rmgte,n,freq,shortTree,mgtepmtable)

r=(1+rmgte)^(1/freq)-1;    
T=n*freq;
dt=1/freq;
mgtepmt=fv*r/(1-1/(1+r)^T);    

optimal_Tree=nan(T+1,T+1);
optimal_Tree(end,1:end)=mgtepmt;
alt_tree=nan(361,361);              
alt_tree(end,1:end)=mgtepmt;
for i=T:-1:1
    for j=i:-1:1
        if i==360
            alt_tree(i,j)=0.5*(alt_tree(i+1,j)+alt_tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(alt_tree(i,j),mgtepmtable(i)) + mgtepmt;   
        elseif i==1
            alt_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(alt_tree(i,j),mgtepmtable(i));
        else
            alt_tree(i,j)=0.5*(optimal_Tree(i+1,j)+optimal_Tree(i+1,j+1))*exp(-shortTree(i,j)*dt);
            optimal_Tree(i,j) = min(alt_tree(i,j),mgtepmtable(i)) + mgtepmt;
        end
    end
end


p0=optimal_Tree(1,1);       
pU=optimal_Tree(2,2);      
pD=optimal_Tree(2,1);    
pUU=optimal_Tree(3,3);
pUD=optimal_Tree(3,2);
pDD=optimal_Tree(3,1);

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
%output=optimal_Tree;