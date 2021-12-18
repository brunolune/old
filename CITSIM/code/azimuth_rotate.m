function [x1,y1,z1] = azimuth_rotate(az,x,y,z)
az=az*2*pi/360
x1=z*sin(az)+x*cos(az);
y1=y;
z1=z*cos(az)-x*sin(az);
end