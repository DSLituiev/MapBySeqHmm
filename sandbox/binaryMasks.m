bN = 4;
diagMask = eye(bN);
b10 = num2str( diagMask*10.^(0:bN-1)');
binMask = uint8(bin2dec(b10'));

bitget(uint8(bin2dec('0101')), 3)

binMask