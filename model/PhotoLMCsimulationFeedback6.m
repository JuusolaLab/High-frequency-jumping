%Script to simulate LMC voltage from photovoltage with feedback

%Parameters for synaptic probability model
LMCIparam.ptau =0.1;% activation time constant
LMCIparam.ntau =7.5;%inacativation time constant
LMCIparam.V0= -55.9483;%Half activation value
LMCIparam.dV =2.5; %Slope of activation


%Length of simulation
time_length=length(Xstim)/samplerate*1000;
%Simulation initial values for pd(0) and pf(0) 
yinitsyn =[0 0]; %Initial values for pd(0) and pf(0) simulation
 time_resol =1/samplerate*1000; %Time resolution in ms
%Parameters for stocastic model
tspan = [0:time_resol:time_length]; % sampling time points
 optionssyn=odeset('RelTol',1e-6,'maxstep',0.5);
 yI =zeros(length(yinitsyn),length(Xstim),6);%pd(t) and pf(t)
 yI(:,1) =yinitsyn;
 SYNInput =zeros(6,length(Xstim));% probability of activate terminal

%Stocastic model

MaxActiveterminals = 390*6;%Number of LMC terminals
 LatencyGamma =[1 2.34e-05];%LMC Latency gamma 
 RefractoryGamma =[0.027	0.585];%LMC Refractory gamma

tt =(0:15)/fs; %Time axis for LMC bump
sstateprob =0.0529; % Tonic release probablility
%Initiliaze vectors
datalength =length(SYNInput);
Activeterminals = zeros(6,datalength);%Number of active terminals
N_Activeterm =zeros(6,1)+MaxActiveterminals/6;%Current active terminals
ReleasedTerminals =zeros(6,datalength);%Number of terminal activated
AfterLatency =zeros(6,datalength);%Number of terminal after latency
RefractorySituation=zeros(6,datalength);%Number of terminal coming back from refractory 
TerminalActivations = zeros(MaxActiveterminals,datalength);%Terminal activations
gLICind =zeros(MaxActiveterminals,datalength+length(tt));%Terminal conductances
gLICm =zeros(1,datalength);%Total histamine conductance
 
%Membrane parameters
paramLMCMembrane.Cm=0.9;%*10^-6 F/cm^2 Specific capacitance
Cm =1.8000e-11;% Membrane capacitance
paramLMCMembrane.Sm = Cm/paramLMCMembrane.Cm/1e-6;%cm^2 Membrane area
 Rm =21.87e6;%ohm Passive resistance
 paramLMCMembrane.gRest= 1/Rm/paramLMCMembrane.Sm;%S/cm^2 
 paramLMCMembrane.gKd =0.006;%S/cm^2 Dr maximal conductance


 paramLMCMembrane.EK =-65;%mV Dr reverse potential
 paramLMCMembrane.VLIC =-85;%mV Histamine channel reverse potetial
 paramLMCMembrane.VRest =-29.3067; %mV resting potential
 %Initial parameters for LMC model
 V =-49.11; %mV starting voltage
 %Kd parameters
 hmidv = -84;% (mV)
 hslope = 8.7 ;%(mV)
 mmidv = -60; %(mV)
 mslope = -5.5;% (mV)
 Kdsift = 20; %(mV) Sift in Kd values

 paramLMCMembrane.Kdsift =Kdsift;
 %Initial Kd parameters
 m0 =1/ ( 1 + exp( (V-(mmidv+Kdsift))/mslope ) );
 h0 =1/ ( 1 + exp( (V-(hmidv+Kdsift))/hslope ) );

%Conductance bump
gBumpparam =[2.16e-13,2.67,4.38e-04]; %LMC conductance bump parameters
[gBump gBumpDuration] = Bumpcalc(gBumpparam(1),gBumpparam(2),gBumpparam(3),tt); %Calculate LMC conductance bump

%Initialise vectors
ILMC =zeros(1,length(gLICm)); %LMC current injection
yinitLMCmemb =[V m0 h0]; %Initial values for LMC membrance simulation
yLMC =zeros(length(yinitLMCmemb),length(gLICm)); % LMC voltage, m,h
yLMC(:,1) =yinitLMCmemb;
 optionsLMCmem=odeset('RelTol',1e-3,'maxstep',0.5);

 %Main loop
for i = 1:datalength
    %Synaptic Input
    for k =1:6
        if(i>1)
            %Simulate differential equation postsynaptic activation
            %probability with feedback
            temp_y =ode15s(@SynapseModel,tspan([(i-1) i]),yI(:,i-1,k),optionssyn,Xstim(k,:),yLMC(1,:),time_resol,LMCIparam,1); 
            
            yI(:,i,k)= temp_y.y(:,end);
            SYNInput(k,i) =sstateprob+(1-sstateprob)*yI(1,i,k);
            if(SYNInput(k,i)>1)
                SYNInput(k,i)=1;

            end

        end
        %Update active postsynaptic sites
        N_Activeterm(k) =N_Activeterm(k)+RefractorySituation(k,i);
        Activeterminals(k,i) =N_Activeterm(k);
        %Calculate number of synaptic activations
        N_ActTerm=binornd(N_Activeterm(k),SYNInput(k,i));
        ReleasedTerminals(k,i) =N_ActTerm;

        N_Activeterm(k) =  N_Activeterm(k)-N_ActTerm;

        if(N_ActTerm >0) %If currently any activations
            %Calculate latency delay
            Latencies = round(gamrnd(LatencyGamma(1),LatencyGamma(2),N_ActTerm,1)*fs);
            Latencies_time = Latencies +i;
            Latencies_time(Latencies_time>datalength)=[];
            if(~isempty(Latencies_time))
                [Nlat,poslat] = histcounts(Latencies_time,'BinMethod','integers');
                AfterLatency(k,poslat(1:(end-1))+0.5) =AfterLatency(k,poslat(1:(end-1))+0.5)+Nlat;
            end
            %Calculate refractory
            RefractroyPeriod = round(gamrnd(RefractoryGamma(1),RefractoryGamma(2),N_ActTerm,1)*fs);

            Refraftory_time = i+RefractroyPeriod;
            Refraftory_time(Refraftory_time>datalength-1)=[];
            if(~isempty(Refraftory_time))
                [Nref,posref] = histcounts(Refraftory_time,'BinMethod','integers');
                 RefractorySituation( k,posref(1:(end-1))+1.5) =RefractorySituation(k, posref(1:(end-1))+1.5)+Nref;
            end
        end
        %Random site activations
        randorder = randperm(MaxActiveterminals/6);
        TerminalActivations((k-1)*MaxActiveterminals/6+randorder(1:AfterLatency(k,i)),i)=1;
    end
  %Calculate histamine conductance 
  gLICind(:,i:(i+length(tt)-1))=gLICind(:,i:(i+length(tt)-1))+TerminalActivations(:,i)*gBump;
  temp_gl =gLICind(:,i);
  temp_gl(temp_gl> max(gBump)) = max(gBump);
   gLICind(:,i)=temp_gl;
  
 gLICm(i) =sum(temp_gl)/paramLMCMembrane.Sm;
    %Simulate LMC membrane model differential equation
    if(i>1)
         temp_y=ode45(@gLICmodel,tspan([(i-1) i]),yLMC(:,i-1),optionsLMCmem,gLICm,ILMC,time_resol,paramLMCMembrane);
         yLMC(:,i)= temp_y.y(:,end);
    end
end

%Synapse differential equation
function dy=SynapseModel(t,y,V,yLMC,time_resol,param,Diameter_ratio)
time_point=round(t/time_resol+1);   % Current stimulus at time t
if time_point>length(V), Vc=0; else Vc=V(time_point); end
if time_point>length(yLMC), YLMC=0;elseif time_point==1, YLMC=yLMC(1); else YLMC=yLMC(time_point-1); end

 
Vc =Vc-7*y(2);%Feedback amount
Diacorr=1;
decay = 1/param.ntau/(1+(0.05/y(1))^2);%Decay back to no synaptic activation
dy =[1/param.ptau/Diacorr*(exp((Vc-param.V0)/param.dV)*(1-y(1)))-decay;%d pd(t)
   1/0.8*(exp(-1/2.5*(YLMC+56))*(1-y(2))-y(2))];%d of(t)

end

%Membrane differential equation
function dy=gLICmodel(t,y,gLICm,Im,time_resol,param)
time_point=round(t/time_resol+1);   % Current stimulus at time t
if time_point>length(Im), I=0; else I=Im(time_point); end
if time_point>length(gLICm), gLIC=0; else gLIC=gLICm(time_point); end
%Calculate Dr channel properties
m =y(2);
h =y(3);
eK = param.EK;% (mV) : Potassium reverse potential
	gbar =  param.gKd;
	
	mmidv = -60; %(mV)
	mslope = -5.5;% (mV)
	mtaumax = 35.6-4.2;% (ms)
    mtaumin =4.2;% (ms)
	mmidvd = -62.6;% (mV)
	msloped = 3.4;% (mV)

	hmidv = -84;% (mV)
	hslope = 8.7 ;%(mV)
	htaumax = 2000 ;%(ms)

	Kdsift = param.Kdsift; %(mV)
    minf = 1/ ( 1 + exp( (y(1)-(mmidv+Kdsift))/mslope ) );
	mtau = mtaumax / (1+ exp( (y(1)-(mmidvd+Kdsift))/msloped))+mtaumin ;
    hinf = 1/ ( 1 + exp( (y(1)-(hmidv+Kdsift))/hslope ) );
	htau = htaumax;
    gKd = gbar*m*h;
    %Differential equations
    dy =[(-y(1)*(param.gRest+gLIC+gKd)/(0.001*param.Cm)+(param.gRest*param.VRest+gLIC*param.VLIC+gKd*eK)/(0.001*param.Cm))+(I/(1000*param.Cm*param.Sm));...
     (minf - m)/mtau;...
     (hinf - h)/htau];


end
