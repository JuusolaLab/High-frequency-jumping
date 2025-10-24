% Contains various photoreceptor model paramemeters
Rm = 4.0239e+07;%Resting resistance
Cm = 1.78e-11;%Capacitance
%Membrane parametrization
param.Cm=0.9;%*10^-6 F/cm^2
param.Sm = Cm/param.Cm/1e-6;%cm^2
param.VRest=-76.8;%Resting reverse potential
param.VK= -85; %mV Potassium reverse potential
param.gShab =10e-3; %S/cm^2 Shab conductance 30
param.gRest= 1/Rm/param.Sm;%S/cm^2
param.VLIC =10;%mV light Nerst

gBumpparam =[2.50e-14,1.24,5.77e-04];%Conductance bump highres 200 Hz bursty BG0


