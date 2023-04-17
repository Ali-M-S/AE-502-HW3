
function [h_eqnx,k_eqnx,p_eqnx,q_eqnx] = ic(a,e,i,M,omega,Omega,w)
mu =1;
deg=180/pi;

%initial conditions
i = deg2rad(i);
M = deg2rad(M);
omega = deg2rad(omega);
Omega = deg2rad(Omega);

n = sqrt(1/(a^3)); % rad/s 

% initial given actions and angles
L0 = sqrt(mu*a);
G0 = L0*sqrt(1 - e^2);
H0 = G0*cos(i);
l0 = M; %n*TU at 0
g0 = omega;
h0 = Omega+omega;

de = [L0 G0 H0 l0 g0 h0 L0 G0 H0 l0 g0 h0];

t0 = 0;
tf = 100;
nout = 10000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset(...
'reltol', 1.e-10, ...
'abstol', 1.e-10);

y0 = de';

[t,y] = ode45(@rates, tspan, y0, options);
%...Assign the time histories mnemonic variable names:
L = y(:,1);
G = y(:,2);
H = y(:,3);
l = y(:,4);
g = y(:,5);
h = y(:,6);

Le = y(:,7);
Ge = y(:,8);
He = y(:,9);
le = y(:,10);
ge = y(:,11);
he = y(:,12);

La = Le + w*(Le.^3).*He;
Ga = Ge;
Ha = He;
la = le - 3*w*((Le.^2).*He.*le);
ga = ge;
ha = he - w*((Le.^3).*le);

h_eqnx = e*sin(ga+ha);
k_eqnx = e*cos(ga+ha);
p_eqnx = tan(i/2)*sin(ha);
q_eqnx = tan(i/2)*cos(ha);



hold on
plot(h_eqnx,k_eqnx,'r','LineWidth',2);
plot(p_eqnx,q_eqnx,'b','LineWidth',2);
legend('h vs k','p vs q', 'Interpreter', 'Latex')
hold off


function dfdt = rates(t,f)

% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

L = f(1);
G = f(2);
H = f(3);
l = f(4);
g = f(5);
h = f(6);

Le = f(7);
Ge = f(8);
He = f(9);
le = f(10);
ge = f(11);
he = f(12);

% vinti pg. 108
Ldot = 0;
Gdot = 0;
Hdot = 0;
ldot = 1/(L^3);       %- w*sqrt(1 - e^2)*cos(i);
gdot = 0; %-w*cos(i);
hdot = w;

Ldot2 = 0;
Gdot2 = 0;
Hdot2 = 0;
ldot2 = 1/(Le^3);       
gdot2 = 0; 
hdot2 = 0;

dfdt = [Ldot Gdot Hdot ldot gdot hdot Ldot2 Gdot2 Hdot2 ldot2 gdot2 hdot2]';


end
end