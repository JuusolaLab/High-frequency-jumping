function [ Data] =Lightmodelsimplewithfeedback2dot(modelparam,Data)
%Lightmodelsimplewithfeedback simulate light absorbtion based on light
%intensity taking account microsaccadic movement

% modelparam model paremeters
% Data input/output

%Set the seed for random generator
rndseed = RandStream.create('mt19937ar','seed', round(cputime*100));
RandStream.setGlobalStream(rndseed);

steps = length(Data.light_series);
n_photoreceptors = size(Data.light_series,2);
if(length(modelparam.N_micro)==1)
    Pro = 1/modelparam.N_micro*ones(modelparam.N_micro,1);
end
%Create outputs
Data.Absorbtions = zeros(steps,n_photoreceptors);% Light absorbtions
Data.AfterLatency =zeros(steps,n_photoreceptors);% Light absorbtions after latence
Data.Refractory_situation =zeros(steps,n_photoreceptors);%Refractory return
Data.OUT =zeros(steps,n_photoreceptors);%Light current
Data.xpos = zeros(steps+1,1);% Rhabdomere position
Data.xposd = zeros(steps+1,1); % Speed of the rhabdomere
% Data.xpos(1) =modelparam.xposstart;
%time series for single bump
fs = modelparam.Fs;
tt =(0:30)/fs;
[BumpShape BumpDuration] = Bumpcalc(modelparam.BumpParam(1),modelparam.BumpParam(2),modelparam.BumpParam(3),tt);


%Main loop
for i = 1:steps
    %Rhabdomere movement calculation
    %Spring constant
    act=sum(Data.AfterLatency(i,:))/modelparam.actdiv;
    
      % spring =modelparam.spring_0+modelparam.spring_coef*(act/modelparam.spring_1_point)^modelparam.spring_n;
    spring =modelparam.spring_0+modelparam.spring_coef*act^modelparam.spring_n/(modelparam.spring_1_point^modelparam.spring_n+act^modelparam.spring_n);
    % dv/dt for Rhabdomere movement
      dvdt = modelparam.Activation_Force_max*act^ modelparam.Activation_Force_n/(modelparam.Activation_Force_1_point^modelparam.Activation_Force_n+act^modelparam.Activation_Force_n)...
         +modelparam.Dampener_coef*(modelparam.Dampener_base^(-1*modelparam.Dampener_exponent*Data.xposd(i))-1)-spring*Data.xpos(i);
    Data.xposd(i+1) =Data.xposd(i)+1000/fs*dvdt;
    % dx/dt for Rhabdomere movement
    Data.xpos(i+1) =Data.xpos(i)+ 1000/fs*Data.xposd(i);
    
    %Calculate new receptive fields
    if(Data.xpos(i)~= Data.xpos(i+1))
        Data.Map = MultipleFields(modelparam.MU+repmat([ Data.xpos(i+1)*modelparam.AngleScale(1) Data.xpos(i+1)*modelparam.AngleScale(2)],size(modelparam.MU,1),1),modelparam.hw-modelparam.hwmd *Data.xpos(i+1),modelparam.amplitude+modelparam.amplitudemd*Data.xpos(i+1),modelparam.xdim,modelparam.ydim,Data.lensposc,Data.lensnormalc,modelparam.mappos, modelparam.mapsize);
    end
    bar_x_pos =modelparam.xstartpos+round((i-1)*modelparam.barspeed); %Dot position
    Barpos =modelparam.barpos;
    Barpos(:,1) =Barpos(:,1)+bar_x_pos; %Dot position

    Data.light_series(i,:)  = BarVideo(Data.Map,Barpos , 1,modelparam.barangle, 1); %Calculate light intensities for dot


    %Loop through ommatidium photoreceptors
    for k =1:n_photoreceptors
        %Photon count absorbed
        N_photon = round(Data.light_series(i,k)*modelparam.LightScale);
        %Calculate number of activate microvilli
        if(length(modelparam.N_micro)>1)
               Active_microvilli = modelparam.N_micro(k)+sum(Data.Refractory_situation(1:i,k)-Data.Absorbtions(1:i,k)); % Microvilli from refractory state to active state
        else
            Active_microvilli = modelparam.N_micro+sum(Data.Refractory_situation(1:i,k)-Data.Absorbtions(1:i,k)); % Microvilli from refractory state to active state
        end
        if(N_photon ~=0 && Active_microvilli ~=0)
             %Claculate number of activated microvilli 
            if(length(modelparam.N_micro)>1)
                   Absorb_all = mnrnd(N_photon,1/modelparam.N_micro(k)*ones(modelparam.N_micro(k),1),1);% all microvilli absorbtions
            else
                Absorb_all = mnrnd(N_photon,Pro,1);% all microvilli absorbtions
            end
            Absorbtions_cur = sum(Absorb_all(1:Active_microvilli)>0);% Absorbtions in active microvilli
            if(Absorbtions_cur >0) %If currently any absorbtions
                Data.Absorbtions(i,k) =Absorbtions_cur;
                %Get the latencies of activated microvilli
                Latencies =  round(gamrnd(modelparam.LatencyDis(1),modelparam.LatencyDis(2),Absorbtions_cur,1)*fs);
                Latencies_time = Latencies +i;
                Latencies_time(Latencies_time>steps)=[];
                if(length(Latencies_time)>0)
                    [Nlat,poslat] = histcounts(Latencies_time,'BinMethod','integers');
                    Data.AfterLatency(poslat(1:(end-1))+0.5,k) =Data.AfterLatency(poslat(1:(end-1))+0.5,k)+Nlat';
                end
            end
        else
            Data.Absorbtions(i,k) = 0;
        end
        if( Data.AfterLatency(i,k)>0) %If start of bump any at the moment
           %      Calculate the refractory of activated microvilli
            RefractroyPeriod = round(gamrnd(modelparam.BumpRefracDis(1),modelparam.BumpRefracDis(2),Data.AfterLatency(i,k),1)*fs);
            Refraftory_time = i+BumpDuration+RefractroyPeriod;
            Refraftory_time(Refraftory_time>steps)=[];
            if(length(Refraftory_time)>0)
                [Nref,posref] = histcounts(Refraftory_time,'BinMethod','integers');
                Data.Refractory_situation( posref(1:(end-1))+0.5,k) =Data.Refractory_situation( posref(1:(end-1))+0.5,k)+Nref';
            end
        end
    end
end
end

