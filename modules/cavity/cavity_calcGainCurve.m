function gain_curve = cavity_calcGainCurve(w,k_a,g_per)
    gain_curve = g_per^2./(g_per^2 + (w-k_a).^2);
end
