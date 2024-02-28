function neg_return = ivtsObjective(ivtsini, VXX, VXZ, VIX, VIX3M, workbook, w)
    [~, neg_return] = ivtsopt(ivtsini, VXX, VXZ, VIX, VIX3M, workbook, w);
end
