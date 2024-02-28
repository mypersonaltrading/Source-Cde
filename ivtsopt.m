function [outputval,neg_return] = ivtsopt(ivtsini,VXX,VXZ,VIX,VIX3M,workbook,w)
%input vix vix3m indices, and vxx vxz(these two could be changed into VIX
%futrure prices etc.
%   if you wish to trade other products, see details below.
VIX.IVTS=VIX.AdjClose./VIX3M.AdjClose;



for i=2:height(VIX)
    if VIX.IVTS(i)<=ivtsini(1);
        VXX.Weight(i)=w(1);%-0.60
        VXZ.Weight(i)=1+VXX.Weight(i);
    elseif VIX.IVTS(i)<=ivtsini(2) & VIX.IVTS(i)>ivtsini(1);
        VXX.Weight(i)=w(2);%-0.32
        VXZ.Weight(i)=1+VXX.Weight(i);
    elseif VIX.IVTS(i)<=ivtsini(3) & VIX.IVTS(i)>ivtsini(2);
        VXX.Weight(i)=w(3);%-0.25
        VXZ.Weight(i)=1+VXX.Weight(i);
    else VIX.IVTS(i)>ivtsini(3);
        VXX.Weight(i)=w(4);%-0.10
        VXZ.Weight(i)=1+VXX.Weight(i);
    end

    %above if loop is used to determine vxx and vxz weight, however, you
    %can change vxx and vxz into VIX 1M future and VIX 3M futures.
    workbook.VXXWadjR(i)=VXX.Weight(i)*VXX.Return(i);
    workbook.VXZWadjR(i)=VXZ.Weight(i)*VXZ.Return(i);
    workbook.momentum(i)=workbook.VXXWadjR(i)+workbook.VXZWadjR(i);
     
end

workbook.CumRet=cumprod(1+workbook.momentum);
outputval=[ivtsini(1),ivtsini(2),ivtsini(3)];
overallreturn=workbook.CumRet(502);
neg_return=-overallreturn;


end