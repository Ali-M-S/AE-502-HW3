clear all; clear;clc

deg=180/pi;
%initial conditions
mu = 1;
a = 1;
e = 0.5;
i = deg2rad(45);
w = 0.01;

n = sqrt(1/(a^3)); % rad/s 

% initial given actions and angles
L0 = sqrt(mu*a);
G0 = L0*sqrt(1 - e^2);
H0 = G0*cos(i);
l0 = 0; %n*TU at 0
g0 = 0;
h0 = 0;

de = [L0 G0 H0 l0 g0 h0];

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


for p = 1:length(l)

% Findnig Eccentric anomaly (initial)
E = l(p,1); %Initial guess for eccentric anomaly E = Mean anomaly
s = 1;
itr = 0; 
while abs(s) > 1e-13  
s = E-e*sin(E) - l(p,1); % rad/s (modified Kepler's equation)
dgdE = 1-e*cos(E);
E_new = E - s/dgdE;
% Update
E = E_new;
itr = itr + 1;
end
E_vec(p,1) = E;
end

% Find True Anomaly from Eccentric Anomaly
 f1 = (2.*atan(sqrt((1 + e) ./ (1 - e)) .* tan(E_vec/2))); 
 kk = E_vec./(2*pi);              %    Case 1       or   Case 2
 k_round = round(kk);   %k_round     = 1               = 0
 f = f1 + k_round*(2*pi);

%find r magnitude from true anomaly
for u=1:length(f)
 r(u)  = a*(1 - e^2)/(1 + e*cos(f(u)));
 
 % X-Y-Z position components
 % from https://www.jstor.org/stable/2635523
 xf(u) = r(u)*cos(f(u)+g(u))*cos(h(u)) - r(u)*sin(f(u)+g(u))*cos(i)*sin(h(u));
 yf(u) = r(u)*cos(f(u)+g(u))*sin(h(u)) + r(u)*sin(f(u)+g(u))*cos(i)*cos(h(u));
 zf(u) = r(u)*sin(f(u)+g(u))*sin(i);

 xg(u) = r(u)*cos(f(u)+g(u))*cos(h(u)) - H(u)/G(u)*r(u)*sin(f(u)+g(u))*sin(h(u));
 yg(u) = r(u)*cos(f(u)+g(u))*sin(h(u)) + H(u)/G(u)*r(u)*sin(f(u)+g(u))*cos(h(u));
 zg(u) = sqrt(G(u)^2 - H(u)^2)/G(u)*r(u)*sin(f(u)+g(u));
end

%inertial frame components (AE 402 lecture 22)
% since angular velocity is constant
xi   = xf.*cos(w.*t') - yf.*sin(w.*t'); %x-component
eta  = xf.*sin(w.*t') + yf.*cos(w.*t'); %y-component
zeta = zf;                              %z-component
figure (1)
hold on
title('Three-Body Trajectory', 'Interpreter', 'Latex')
xlabel('x-position (AU)')
ylabel('y-position (AU)')
zlabel('z-position (AU)')
axis auto
grid on
view(-0.5,0);
plot3(xf, yf, zf, 'r','LineWidth',1)
legend('Rotational frame', 'Interpreter', 'Latex')
hold off

figure(2)
hold on
title('Three-Body Trajectory', 'Interpreter', 'Latex')
xlabel('x-position (AU)')
ylabel('y-position (AU)')
zlabel('z-position (AU)')
axis auto
grid on
view(-0.5,0);
% Plotting Trajectory
plot3(xf, yf, zf, 'r','LineWidth',1)
plot3(xi, eta, zeta, 'k','LineWidth',1)
legend('Rotational frame','Inertial frame', 'Interpreter', 'Latex')
hold off

function dfdt = rates(t,f)

% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mu = 1;
a = 1;
e = 0.5;
i = deg2rad(45);
n = sqrt(1/(a^3)); % rad/s 
w = 0.01;

L = f(1);
G = f(2);
H = f(3);
l = f(4);
g = f(5);
h = f(6);

% vinti pg. 108
Ldot = 0;
Gdot = 0;
Hdot = 0;
ldot = 1/(L^3);       %- w*sqrt(1 - e^2)*cos(i);
gdot = 0; %-w*cos(i);
hdot = w;

dfdt = [Ldot Gdot Hdot ldot gdot hdot]';

end


