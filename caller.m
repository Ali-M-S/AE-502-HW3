clear all; clear;clc

%Please choose between "e"   for eccentricity
%                    "Omega" for Ω
%                      "w"   for rotation frequency


%%%%-------------
conditions = "e";
%%%%-------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if conditions == "e"
for e=0:0.1:0.9
%                                  (a, e, i,M,omega,Omega,w)
[h_eqnx,k_eqnx,p_eqnx,q_eqnx] = ic(1,e,45,0,0,0,0.01);
hold on
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif conditions == "Omega"
for Omega=[10 45 90 80] %degrees
%                                  (a, e, i,M,omega,Omega,w)
[h_eqnx,k_eqnx,p_eqnx,q_eqnx] = ic(1,0.5,45,0,0,Omega,0.01);
hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif conditions == "w"
for w=[0.02, 0.1, 0.5]
%                                  (a, e, i,M,omega,Omega,w)
[h_eqnx,k_eqnx,p_eqnx,q_eqnx] = ic(1,0.5,45,0,0,0,w);
hold on
end
end
