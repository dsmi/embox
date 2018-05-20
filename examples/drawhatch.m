function BH = drawhatch(B0, hp, hw, ha, ho)

n = round(1.2/hw);
l = 1.2;
BH = B0;

for i = -n : n

  % horz
  xc =            hp*i * sin(ha);
  yc = 0.5 + ho + hp*i * cos(ha);
  x0 = xc + l * cos(ha);
  y0 = yc - l * sin(ha);
  x1 = xc - l * cos(ha);
  y1 = yc + l * sin(ha);
  BH = linefromto(BH, x0, y0, x1, y1, hw);

  % and vert
  xc =            hp*i * cos(ha);
  yc = 0.5 + ho - hp*i * sin(ha);
  x0 = xc + l * sin(ha);
  y0 = yc + l * cos(ha);
  x1 = xc - l * sin(ha);
  y1 = yc - l * cos(ha);
  BH = linefromto(BH, x0, y0, x1, y1, hw);

end
