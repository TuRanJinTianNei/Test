function [classes, counts] = identify_sig(y, centers)
    classes = zeros(1, length(y));
    counts = zeros(1, length(centers));
    for i = 1:length(y)
        tmp = abs(centers - y(i))./centers.^2;
        [~, classes(i)] = min(tmp);
        counts(classes(i)) = counts(classes(i)) + 1;
    end
end