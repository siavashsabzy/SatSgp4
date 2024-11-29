clear all
clc
%%
line1 = "1 24876U 97035A   24310.87505262  .00000074  00000+0  00000+0 0  9995 ";
line2 = "2 24876  55.7186 121.5906 0085632  54.8461 305.9738  2.00561646200150 ";
name = "NAVSTAR 43 (USA 132)    ";
satrec = SatSgp4(char(line1), char(line2), char(name));
eciPosVel_meter = satrec.propagate(2460644.85416667); % externally validated