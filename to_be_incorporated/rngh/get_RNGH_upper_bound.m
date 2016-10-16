function RNGH_u = get_RNGH_upper_bound(lambda,gpar,kappa)
    R = (lambda.^2 - 8*lambda - 6*gpar*lambda + gpar^2).^0.5;   %R/gper
    RNGH_u = (0.5*gpar*(3*lambda - gpar + R)).^0.5.*(1-2*kappa./(lambda - 2 - gpar + R));
end