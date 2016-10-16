function psteps = pump_getPumpSteps(pump,pump0,pump1)   %all quantities are *NOT* normalized by threshold. ie this is 'pump' not 'beta'
    [~,ind0] = helpers_getClosestMatch(pump0,pump);
    [~,ind1] = helpers_getClosestMatch(pump1,pump);
    
    psteps = ind0:ind1;
end