function [match,index] = helpers_getClosestMatch(seed,matches)
    [~,index] = min(abs(matches - seed));
    match = matches(index);
end