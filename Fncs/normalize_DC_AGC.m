function y = normalize_DC_AGC(x,kDC,kAGC)
% kDC  = 1e4;                                                     % memory for DC removal
% kAGC = 1e4;                                                     % memory for Automatic Gain Control

% Normalize DC and AGC
x = x-movmean(x,kDC,1);                                                    % remove DC
y = x./movstd(x,kAGC,1);                                                   % 归一化

end