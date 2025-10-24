%Simulates section of bursty 200Hz stimulus for Figure 4
%@Jouni Takalo
load Bursty_stim_photo_hres_200.mat;% load for LMC 200 Hz bursty
ImpulseInput = ImpulseInput(1:8000); %Simulate 2 first repeates of stimulus
MultipleRhabdomerewithScreen %Simulate the photoreceptor quantal sampling process with ray tracing and microsaccades 
CombinedAfterLatency %Simulates photoreseptor voltage 
Xstim =Photovoltage;
PhotoLMCsimulationFeedback6%LMC simulations for LMC voltage with feedback 
figure
subplot(3,1,1)
plot(tspan(4002:end),ImpulseInput(4001:end))%Plot stimulus
title('Light Intensity')
ylabel('arbitary light intensity')
subplot(3,1,2)

plot(tspan(4002:end),Photovoltage(:,4001:end))%Plot six photoreceptor voltages for second repeat
title('Photoreceptor voltage')
ylabel('Membrane voltage (mV)')
xlabel('time (ms)')
subplot(3,1,3)
plot(tspan(4002:end),yLMC(1,4001:end))%Plot LMC voltage for second repeat
title('LMC voltage')
ylabel('Membrane voltage (mV)')
xlabel('time (ms)')
