function xdot = myOde(t,x)

k = 398600;

xdot(1) = x(2); 
xdot(2) = x(1)*x(4)^2 - k/(x(1)^2); 
xdot(3) = x(4); 
xdot(4) = -2*x(4)*x(2)/x(1); 

xdot = [xdot(1); xdot(2); xdot(3); xdot(4)];