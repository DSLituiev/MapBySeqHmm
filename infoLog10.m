function y = infoLog10(logX)

logX = logX - calcMarginal(logX);

 y = -sum( (log2(10)*logX).*10.^logX, 1 );

