function Mdata =Photomembrane(AfterLatency,gBumpparam,param,samplerate)
%Simulate photoreceptor membrane voltage from absorbed light
%AfterLatency ligth activations after latensy
%gBumpparam conductance bump parameters
%param cell membrane parameters
%samplerate simulation sample rate

tt =(0:30)/samplerate;%Conductance bump time axis
%Calculate conductance bump
[gBump gBumpDuration] = Bumpcalc(gBumpparam(1),gBumpparam(2),gBumpparam(3),tt);
gLICm =conv(AfterLatency,gBump)/param.Sm;%Calculate light conductance
I =zeros(1,length(gLICm));% Current imput for membrane
time_resol =1/samplerate*1000;%time resolution in ms
time_length=length(gLICm)/samplerate*1000;%Length of simulation in ms
%Shab ion channel
h0=1/(1+exp((-25.7-param.VRest)/-6.4));     % start point at Shab inact. curve
n0=(1/(1+exp((-1-param.VRest)/9.1)))^(1/2); % start point at Shab act. curve
%Initialise parameters
yinit =[param.VRest h0 n0];
tspan = [0:time_resol:time_length];   
 options=odeset('RelTol',1e-4,'maxstep',0.5);
 %Simulate differential equation for membrane voltage
 [t,y]=ode15s(@gLICmodel,tspan,yinit,options,gLICm,I,time_resol,param);
 %Make output
Mdata.t =t;%time
Mdata.y =y;%Voltage and Shab state response
Mdata.gLICm =gLICm;%Light conductance
Mdata.I =I;%added current

end
%Membrane differential equation
function dy=gLICmodel(t,y,gLICm,Im,time_resol,param)
time_point=round(t/time_resol+1);   % Current stimulus at time t
if time_point>length(Im), I=0; else I=Im(time_point); end
if time_point>length(gLICm), gLIC=0; else gLIC=gLICm(time_point); end

VK=param.VK; %mV Potassium 
gKsm=param.gShab; %Shab conductance
%Shab 
corr =1;
% ------- Inactivation: V50=-25.7, s=-6.4; tau=1200 ms
ah=(1/(1+exp((-25.7-(y(1)))/-6.4)))/(corr*1200/1.35);
bh=(1-(1/(1+exp((-25.7-(y(1)))/-6.4))))/(corr*1200/1.35);
% ------- Activation: V50=-1; s=9.1;

an=((1/(1+exp((-1-(y(1)))/9.1)))^(1/2))/(corr/(1.35*(0.116258*exp((-y(1)...
    -25.6551)/32.1933)+0.00659219*(-y(1)-23.8032)/(exp((-y(1)-23.8032)/...
    1.34548)-1))));
bn=((1-(1/(1+exp((-1-(y(1)))/9.1)))^(1/2)))/(corr/(1.35*(0.116258*exp((-y(1)...
    -25.6551)/32.1933)+0.00659219*(-y(1)-23.8032)/(exp((-y(1)-23.8032)/...
    1.34548)-1))));
gKs=(y(3)^2)*y(2)*gKsm; % Shab (or delayed rectifier)

dy =[(-y(1)*(param.gRest+gLIC+gKs)/(0.001*param.Cm)+(param.gRest*param.VRest+gLIC*param.VLIC+gKs*VK)/(0.001*param.Cm))+(I/(1000*param.Cm*param.Sm)); ...
     -y(2)*(ah+bh)+ah;...
    -y(3)*(an+bn)+an];
end




 
