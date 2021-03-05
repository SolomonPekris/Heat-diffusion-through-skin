clear
T = SkinTempDiffusion1D(0.01,"no",3,393.15,"no");
time = 0:50/0.01;
plot(time*0.01,T(12+1,:))
hold on
plot(time*0.01,T(18+1,:))
plot(time*0.01,T(24+1,:))

T = SkinTempDiffusion1D(0.01,"yes",3,393.15,"no");
time = 0:50/0.01;
plot(time*0.01,T(12+1,:))
hold on
plot(time*0.01,T(18+1,:))
plot(time*0.01,T(24+1,:))

xlabel('Time (s)')
ylabel('Temperature (K)')
legend ('x = 0.00333','x = 0.005','x = 0.00667')

hold off