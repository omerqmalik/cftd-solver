function S_structdata = structdata_load(cav_dir,num,data_id,data_type,chunk,varargin)
    S_setupdata = setup_loadParameters(cav_dir,num);
    calc_times  = benchmark_loadTimes([cav_dir '/' S_setupdata.times_dir]);
    
    if isempty(varargin)
        psteps = 1:length(S_setupdata.pump);
    else
        pump0 = S_setupdata.th*varargin{1};
        if length(varargin) == 1
            pump1 = pump0;
        elseif length(varargin) == 2
            pump1 = S_setupdata.th*varargin{2};
        end
        psteps = pump_getPumpSteps(S_setupdata.pump,pump0,pump1);
    end
    psteps = psteps(sum(calc_times(:,psteps),1) > 0);
    
    [t,Y]  = userdata_load([cav_dir '/' S_setupdata.data_dir],data_id,data_type,psteps);
    if ~strcmp(data_id,'D')
        t = helpers_truncVector(t.',chunk);
        Y = helpers_truncVector(Y,chunk);
    end

    S_structdata   = structdata_getGeneralProperties(data_id,data_type,S_setupdata,calc_times,psteps);
    S_structdata.t = t;
    S_structdata.Y = Y;
    S_structdata.psteps = psteps;
end