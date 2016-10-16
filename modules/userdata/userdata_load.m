function [t,Y] = userdata_load(data_dir,data_id,data_type,psteps)
    fprintf('Loading p = %g of %g...',1,length(psteps));
    load(rawdata_getFileName(data_dir,data_id,data_type,1));
    load(rawdata_getFileName(data_dir,data_id,data_type,'t'));
    fprintf('done.\n');

    if length(psteps) == 1
        Y = Yint;
    else
        Y = zeros(length(t),size(Yint,2),length(psteps));
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
end