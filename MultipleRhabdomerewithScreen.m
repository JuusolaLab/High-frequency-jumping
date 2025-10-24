%Simulate photorecptor light absorbtion on 6 ommatidium
%@Jouni Takalo

%centre lens position and normal
Diameters = [1.5 1.3 1.5 1.3 1.3 1.5 1.1];% (um) Rhabdomere diameters

lenspos =[0 0 0;0 0 0 ;0 0 0; 0 0 0; 0 0 0;0 0 0;0 0 0];% (um) Lens positions
lensnormal = [0 0 1; 0 0 1; 0 0 1; 0 0 1;0 0 1; 0 0 1;0 0 1];%Lens normals

ommatidiumpairs =[1,0;1,-1;1,-2;0,-1;-1,0;-1,1;0,0];%R1-R7/R8 pairings

%parameters of the lens array
lensdist = 22.2;% (um) Lens diameter
lensangle= 3.1/180*pi;% (rad) Interommatidial angle
lensanglexx=lensangle*sqrt(3)/2;
lensanglexy=lensangle/2;
lensnormal =lensnormal./repmat(vecnorm(lensnormal,2,2),1,3);
lensnormbefore = lensnormal;
lensradius = lensdist/lensangle; %radius of eye

%Rhabdomere parameters in rest
rotangle = atan2(0.1*6,1.48);
RotM =[cos(rotangle) -sin(rotangle); sin(rotangle) cos(rotangle)];
degperum =-1.86; %(deg/um position relation to angle)
MU = RotM*[1.19 -1.86; 1.47 -0.19; 1.65 1.44; 0.1 1.48; -1.37 1.73; -1.61 0.13; 0 0]'*degperum;% (um) Rhabdomere angles
modelparam.MU =MU';

modelparam.AngleScale=[-1.86/sqrt(2) -1.86/sqrt(2)];%receptive field micrsaccade movements
modelparam.hw = [2.62 2.10 2.62 2.10 2.13 2.57 2.34];% receptive field half widths
modelparam.hwmd =[0.17 0.16 0.17 0.16 0.16 0.17 0.33];%receptive field halfwidths relation to microsaccades 
modelparam.amplitude = [41.7; 38.9; 41.7;38.9;35.5;43.2;33.0 ];%Maximal absorbtion amplitude
modelparam.amplitudemd =[2.7 3.7 2.7 3.7 3.7 2.7 7.3];%Maximal absorbtion amplitude relation to microsaccades


%Virtual screen parameters
%Rays casted to Screen
modelparam.xdim =-10:0.5:10; %(deg)
modelparam.ydim =-10:0.5:10;%(deg)

%3 deg field
modelparam.distmult =0; %Multipier moving closer and farther
%(um) Screen position 
modelparam.mappos =[-8000/2 8000/2 152000 ;8000/2 8000/2 152000; -8000/2 -8000/2 152000]/2^modelparam.distmult; 
modelparam.mapsize =[9 9];%Number of pixels in screen

modelparam.N_micro = round(54000*Diameters/mean(Diameters(1:6)));%Number of microvilli
modelparam.Fs =2000;%Samplerate
samprate =modelparam.Fs;

%Movement model
%Activation model
modelparam.Activation_Force_n=1;%Activation force n
modelparam.Activation_Force_max =9.8068e-04; %Maximal activation force
modelparam.Activation_Force_1_point =15.3; %Activation half point

modelparam.actdiv =7; %Single or 7 photoreceptors activate
xposinitial =0;%(um) Initial microsaccade position
%Dampener force
modelparam.Dampener_coef =0.00024;
modelparam.Dampener_base = 2;
modelparam.Dampener_exponent = 2100;
% %Spring force
modelparam.spring_0 =0.0001; % without activation spring constant
modelparam.spring_coef =0.001;%spring constant activation coeffisiant
modelparam.spring_1_point =15.3;%half activation value for spring constant
modelparam.spring_n =1;%Spring constant multiplier

%Stocastic parameters in case of taking original photoreceptor
%Latency parameters
 modelparam.LatencyDis = [10.4 0.00037]; 
%Refraction parameters
modelparam.BumpRefracDis =[10.5,0.0034];

%Voltage bump params
modelparam.BumpParam =[0.0011,2.367,7.08e-04];

%Final data stucture
Data = [];

for i = 1:6
    %Ommatidium position
    dlensy =ommatidiumpairs(i,2);
    dlensx =ommatidiumpairs(i,1);
    index =num2str(i);
    %Name in cell
    ci = ['l' index];
    %Lens array position
    Data.(ci).dlensy =ommatidiumpairs(i,2);
    Data.(ci).dlensx =ommatidiumpairs(i,1);   
    

    %Bursty timeseries taking account Poisson nature of light
    
    Data.(ci).light_series = poissrnd(repmat( ImpulseInput(1:8000),1,7))*0.92;
    datalength =length(Data.(ci).light_series);
    Data.(ci).Mapvalue = zeros(length(Data.(ci).light_series),7);
    %Various paramenters for lens

    Data.(ci).Rx = [cos(lensanglexx*dlensx) 0 sin(lensanglexx*dlensx); 0 1 0; -sin(lensanglexx*dlensx) 0 cos(lensanglexx*dlensx)];
    Data.(ci).Ry =[1 0 0; 0 cos(lensangle*dlensy+lensanglexy*dlensx) -sin(lensangle*dlensy+lensanglexy*dlensx); 0 sin(lensangle*dlensy+lensanglexy*dlensx) cos(lensangle*dlensy+lensanglexy*dlensx)];
    Data.(ci).lensnormalc =  Data.(ci).Ry* Data.(ci).Rx*lensnormal';
    Data.(ci).lensnormalc = Data.(ci).lensnormalc';
    Data.(ci).lensposc = lenspos+lensradius*( Data.(ci).lensnormalc-lensnormbefore);
    %Calculate initial values for receptive fields
   
    [ Data.(ci).Map] = MultipleFields(modelparam.MU,modelparam.hw,modelparam.amplitude,modelparam.xdim,modelparam.ydim,Data.(ci).lensposc, Data.(ci).lensnormalc,modelparam.mappos, modelparam.mapsize);
end
%Calculate normalisation factor between number of photons and receptive
%fields
LightScale =0;
for i = 1:6
     index =num2str(i);
    %Name in cell
    ci = ['l' index];
    LightScale = LightScale +sum(Data.(ci).Map(:,:,i),[1 2]);

end

modelparam.LightScale =6/LightScale;

AfterLatency_a =zeros(6,datalength);%Number of quantal samples after latency
Datafields =fieldnames(Data);
Data = struct2cell(Data); 
parfor k =1:6
     index =num2str(k);
         ci = ['l' index];
        Data{k} =Lightmodelsimplewithfeedback2(modelparam,Data{k});% Simulate
         AfterLatency_a(k,:) =Data{k}.AfterLatency(1:datalength,k);
%     end
end

Data = cell2struct(Data,Datafields,1);
