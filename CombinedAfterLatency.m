%Script for calculating LMC response after photoreceptor latency
PhotoParameters;% Photoreceptor parameters
N_microvilli = 54000;%Number of microvilli
Diameters = [1.5 1.3 1.5 1.3 1.3 1.5 1.1];%Rhabdomere diameters
fs =2000;
samplerate =fs;
Photovoltage =zeros(size(AfterLatency_a,1),datalength);%R1-R6 photoreceptor voltages
tt =(0:30)/fs;% time for conductance quantal response
%Setup the loop
gShabOrig=param.gShab;
mDiameters =mean(Diameters(1:6));
Datafields =fieldnames(Data);
Data = struct2cell(Data); 
tparam =param;
clear param;
for k =1:6
    %Calcluate photoreceptor voltage
    AfterLatency =AfterLatency_a(k,:)';
    %Because of different photoreceptor diameters adjust parameters
    Diameterratio =Diameters(k)/mDiameters;
      param =tparam;
    param.Cm=0.9;%*10^-6 F/cm^2
      param.Sm = Cm*(0.25+0.75*(Diameterratio^2))/param.Cm/1e-6;%cm^2
      param.gRest= 1/Rm/Cm*param.Cm*1e-6/(0.25+0.75*(Diameterratio^2));%S/cm^2 %Ratio between body and rhabdomere areas (1:3) also expecting that the ratio between conductance density (3:1)
      param.gShab=gShabOrig/(0.25+0.75*(Diameterratio^2)); %Ratio between body and rhabdomere areas (1:3)

     Mdata =Photomembrane(AfterLatency,gBumpparam,param,samplerate);% Simulate photoreceptor voltage
    
    Data{k}.voltage =Mdata.y(1:datalength,1);
    Photovoltage(k,:) =Mdata.y(1:datalength,1);
end
Data = cell2struct(Data,Datafields,1);

