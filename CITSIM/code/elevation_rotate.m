function [x3,y3,z3] = elevation_rotate(el,x,y,z)
figure
hold on
quiver3(0,0,0,x,y,z,'k-','LineWidth',1,'MaxHeadSize',0.3)
az=-atan2(z,x)
el=el*2*pi/360;
x1=x*cos(az)-z*sin(az)
y1=y;
z1=x*sin(az)+z*cos(az)
quiver3(0,0,0,x1,y1,z1,'b-','LineWidth',1,'MaxHeadSize',0.3)
x2=x1*cos(el)-y1*sin(el)
y2=x1*sin(el)+y1*cos(el)
z2=z1
quiver3(0,0,0,x2,y2,z2,'g-','LineWidth',1,'MaxHeadSize',0.3)
x3=x2*cos(az)+z2*sin(az)
y3=y2
z3=-x2*sin(az)+z2*cos(az)
quiver3(0,0,0,x3,y3,z3,'r-','LineWidth',1,'MaxHeadSize',0.3)
xlabel('x')
ylabel('y')
zlabel('z')
end