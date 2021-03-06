function [t,Y] = userdata_load(data_dir,data_id,data_type,psteps)
    fprintf(['Loading data for ' data_id '_' data_type '\n']);
    fprintf('Loading p = %g of %g...',1,length(psteps));
    load(rawdata_getFileName(data_dir,data_id,data_type,psteps(1)));
    load(rawdata_getFileName(data_dir,data_id,data_type,'t'));
    fprintf('done.\n');

    if length(psteps) == 1
        Y = Yint;
    else
        Y = zeros(size(Yint,1),size(Yint,2),length(psteps));
        Y(:,:,1) = Yint;

        index = 2;
        for i = psteps(2:end)
            fprintf('Loading p = %g of %g...',index,length(psteps));
            load(rawdata_getFileName(data_dir,data_id,data_type,i),'Yint');
            fprintf('done.\n');
            Y(:,:,index) = Yint;
            index = index + 1;
        end
    end
    fprintf('\n');
    
    if strcmp(data_type,'field')
        Y = squeeze(Y);
    end
end