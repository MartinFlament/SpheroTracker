function ip = stretchHistogram(ip, saturated)
    stats = imhist(ip);
    numPixels = numel(ip);
    threshold = numPixels * saturated / 200
    
    % Compute cumulative histogram
    cumHist = cumsum(stats);
    
    % Find min and max based on threshold
    hmin = find(cumHist > threshold, 1, 'first');
    hmax = find(cumHist < (numPixels - threshold), 1, 'last');

    % Rescale
    if hmax > hmin
        ip = imadjust(ip, [hmin/255 hmax/255], []);
    end
end