function setup_dispOverlapIntegralIndex(indexmat,i,j)
    if nargin < 2   %display whole matrix of indices
        s = [];
        N = size(indexmat,1);
        for i = 1:N
            for j = 1:N
                if j < N
                    s = [s strcat(num2str(indexmat(i,j,1),'%2d'),'-',num2str(indexmat(i,j,2),'%2d'),'-',num2str(indexmat(i,j,3),'%2d'),'-',num2str(indexmat(i,j,4),'%2d'),'\t\t')];
                else
                    s = [s strcat(num2str(indexmat(i,j,1),'%2d'),'-',num2str(indexmat(i,j,2),'%2d'),'-',num2str(indexmat(i,j,3),'%2d'),'-',num2str(indexmat(i,j,4),'%2d'),'\n')];
                end
            end
        end
    else    %display just one index
        s = [];
        for r = i
            for c = j
                if c < max(j)
                    s = [s strcat(num2str(indexmat(r,c,1),'%2d'),'-',num2str(indexmat(r,c,2),'%2d'),'-',num2str(indexmat(r,c,3),'%2d'),'-',num2str(indexmat(r,c,4),'%2d'),'\t\t')];
                else
                    s = [s strcat(num2str(indexmat(r,c,1),'%2d'),'-',num2str(indexmat(r,c,2),'%2d'),'-',num2str(indexmat(r,c,3),'%2d'),'-',num2str(indexmat(r,c,4),'%2d'),'\n')];
                end
            end
        end
    end
    fprintf(s);
end