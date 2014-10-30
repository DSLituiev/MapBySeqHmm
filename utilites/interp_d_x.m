function Pos_cM_x = interp_d_x(ChrData, chr,x)

CurrChrData = ChrData{1,chr};

Pos_cM_x = 0.01.*(interp1(CurrChrData(:,1),CurrChrData(:,2),double(x),'cubic','extrap'));

