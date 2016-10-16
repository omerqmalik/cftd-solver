function Yint = rawdata_interpY(T,Y,t)
    leN = size(Y,2);
    Yint = zeros(length(t),leN);
    for i = 1:leN
        [T_this,Y_this] = rawdata_removeDuplicates(T,Y(:,i));
        Yint(:,i) = rawdata_interpComplex(T_this,Y_this,t);
    end
end