function sol = plot_RNGH_window_optim(gper,kappa,lambda)
    lambda_guess = 8;
    alpha = 1.904;
    
    figure;
    col = colormap(autumn(length(gper)));
    
    if length(kappa) == 1
        kappa = ones(1,length(gper))*kappa;
    end
    kappa = kappa./gper;
    
    sol = zeros(length(gper),3);
    for i = 1:length(gper)
        sol(i,1:2) = get_RNGH_optim(kappa(i),alpha/gper(i),lambda_guess,1);
        %check solution
        sol(i,3) = alpha/gper(i) - get_RNGH_upper_bound(sol(i,1),sol(i,2),kappa(i));

        gpar = sol(i,2);
        [alpha_min,alpha_max] = get_RNGH_window(gpar,kappa(i),lambda);
        plot(lambda,real(alpha_min),'color',col(i,:)); hold on;
        plot(lambda,real(alpha_max),'color',col(i,:));
        plot(lambda,imag(alpha_max),'color',col(i,:));
        plot(lambda,imag(alpha_min),'color',col(i,:));
    end
    
    yL = get(gca,'ylim');
    for i = 1:length(gper)
        line([sol(i,1) sol(i,1)],yL,'color',col(i,:));
        line(get(gca,'xlim'),alpha/gper(i)*[1 1]);
    end
end

function [alpha_min,alpha_max] = get_RNGH_window(gpar,kappa,lambda)
    alpha_min = get_RNGH_lower_bound(lambda,gpar,kappa);
    alpha_max = get_RNGH_upper_bound(lambda,gpar,kappa);
end