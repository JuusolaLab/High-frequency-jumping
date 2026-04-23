%Simulates section of double dot stimulus for Figure 6
%@Jouni Takalo
MultipleRhabdomerewithScreenDot %Simulate the photoreceptor quantal sampling process with ray tracing and microsaccades 
CombinedAfterLatency %Simulates photoreseptor voltage 
Xstim =Photovoltage;
PhotoLMCsimulationFeedback6%LMC simulations for LMC voltage with feedback 
figure
subplot(3,1,1)
plot(tspan(1:end-1),LightSeries)%Plot stimulus
title('Light Intensity')
ylabel('arbitary light intensity')
subplot(3,1,2)

plot(tspan(1:end-1),Photovoltage)%Plot six photoreceptor voltages for second repeat
title('Photorec eptor voltage')
ylabel('Membranevoltage (mV)')
xlabel('time (ms)')
subplot(3,1,3)
plot(tspan(1:end-1),yLMC(1,:))%Plot LMC voltage for second repeat
title('LMC voltage')
ylabel('Membrane voltage (mV)')
xlabel('time (ms)')
