function [foundSol,thSol,tvesSol,tvas_prevk] = salt_findCrossoverSol(tvas_thisk,tves_thisk,tvas_prevk)

foundSol = 0;
thSol = 0;
tvesSol = 0;

for j=1:length(tvas_prevk)     %For each eigenvalue, because any one of the eigenvalues could have crossed over
    [~,c]=min(abs(tvas_thisk-tvas_prevk(j)));         %The closest (c) eigenvalue of previous iteration is the same eigenvalue as this iteration
    if (imag(tvas_thisk(c))*imag(tvas_prevk(j))<0)    %A sign change of 'flowing' eigenvalues indicates proximity to the real axis
        foundSol = 1;                                 %Real axis crossed -> solution found
        thSol = 1/real(tvas_thisk(c));
        tvesSol = tves_thisk(:,c);
    end
    tvas_prevk(j) = tvas_thisk(c);            %Previous now becomes current and the order is preserved
end