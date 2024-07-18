%{
Author: Chinmaya Joshi
Organisation: NED3 - University of Arkansas
ASTM D5470 Thermal Resistance Tester
Steady-State Analysis
distance between each point for distances = 1E-06
%}

clc; clear; close all
%Use 'clear' and 'clc' commands after every run

%Input Test Data Filename below:
TestDataFilename = 'TRT-222'; %Test data filenames are saved as TRT-0XX
%Input Sample Thickness below:
ySample = 0.00043   ; %Thickness of sample in m
SSTime = 1400000; %SSTime is the point where the smoothened curve becomes relatively flat i.e. reaches steady state
ySampleLine = linspace(0,ySample, ySample*10^6);

%Analysis Code Starts:
TestNumber = TestDataFilename(1,(5:7));
FileExtension = '.txt';
DataInput = append(TestDataFilename,FileExtension);

TRT = load(DataInput); 
% TRT = readmatrix('Test-100_example.txt','NumHeaderLines',23);
% 'TRT-007.txt','NumHeaderLines',23);
%pulling data from MATLAB workspace file
Data = TRT;
t = Data(:,1); %Time of the experiment in s with 1000Hz sample rate

%Cold side temperatures
C1 = Data(:,7);
C2 = Data(:,6);
C3 = Data(:,5);

%Hot side temperatures
H1 = Data(:,2);
H2 = Data(:,3);
H3 = Data(:,4);

%Smoothening the noise in the data
sH1 = smoothdata(H1);
sH2 = smoothdata(H2);
sH3 = smoothdata(H3);

sC1 = smoothdata(C1);
sC2 = smoothdata(C2);
sC3 = smoothdata(C3);


%Steady state temperatures
sstH1 = mean(sH1(SSTime,:));
sstH2 = mean(sH2(SSTime,:));
sstH3 = mean(sH3(SSTime,:));

sstC1 = mean(sC1(SSTime,:));
sstC2 = mean(sC2(SSTime,:));
sstC3 = mean(sC3(SSTime,:));

%Required test values
ySurface = 0.0044;
yThermocouple = 0.0136;
ySpacing = [0.0044, 0.0180, 0.0316];
n = 3; %number of thermocouplees/datapoints on each side, required for least squares regression

areaSample = 0.016*0.016; %Area of sample in m^2
kAl = 167; %6061 Aluminium Thermal Conductivity

%Calculating temperature gradients between thermocouples

HotSST = [
sstH1 0.0044
sstH2 0.0180
sstH3 0.0316
]; 

ColdSST = [
sstC1 0.0044
sstC2 0.0180
sstC3 0.0316

];

HotFit = LinReg(HotSST);
HotData =plotData(HotFit, HotSST);

ColdFit = LinReg(ColdSST);
ColdData =plotData(ColdFit, ColdSST);

HTslope = HotFit(1,1);
CTslope = ColdFit(1,1);

HTintercept = HotFit(2,1);
CTintercept = ColdFit(2,1);

%{
HTgrad12 = (sstH2 - sstH1)/yThermocouple;
HTgrad23 = (sstH3 - sstH2)/yThermocouple;


HTgrad13 = (sstH3 - sstH1)/(2*yThermocouple); % Avg Temperature Gradient Hot Side

CTgrad32 = (sstC2 - sstC3)/yThermocouple;
CTgrad21 = (sstC1 - sstC2)/yThermocouple;


CTgrad31 = (sstC1 - sstC3)/(2*yThermocouple); %Avg Temperature Gradient Cold Side
%}

%Calculating Heat Flux on both sides

HotFlux = -kAl*HTslope;
ColdFlux = -kAl*HTslope;
AvgFlux = (HotFlux + ColdFlux)/2;
%AvgFlux = HotFlux;
%{
%Least Squares Regression
SST = [sstH1, sstH2, sstH3; sstC1, sstC2, sstC3]; %array of Steady State Temperatures
HT = SST(1,:); %Hot Side Steady State Temps
CT = SST(2,:); %Cold Side Steady State temps

%finding the slopes and intercepts for LSR
sumy = sum(ySpacing);
sumy2 = sum(ySpacing.*ySpacing);
sumHT = sum(HT); 
sumCT = sum(CT);

sumyHT = sum(ySpacing.*HT);
sumyCT = sum(ySpacing.*CT);

%Slopes and Intercepts
Hslope = (n*sumyHT - sumHT*sumy)/(n*sumy2 - sumy^2);
Cslope = (n*sumyCT - sumCT*sumy)/(n*sumy2 - sumy^2);

Hintercept = (sumHT - Hslope*sumy)/n;
Cintercept = (sumCT - Cslope*sumy)/n;
%}
%Using slopes and intercepts from LSR to extrapolate temps on test surfaces
length = linspace(0,0.0360,36000); %length of 1 aluminium block
HTL = HTslope.*length + HTintercept;
CTL = CTslope.*length + CTintercept;

reversedHTL = fliplr(HTL);

slopeSample = (HTL(1,36000) - CTL(1,36000))/ySample;
sampleTemps = slopeSample.*ySampleLine + CTL(1,36000);

totalLength = (2 * 0.0360) + ySample;
experimentalLength = linspace(0,totalLength,totalLength*10^6);
setupPoints = 36000*2 + totalLength*10^6;
setupLength = setupPoints/10^6;

HTLsize = size(HTL,2); 
CTLsize = size(CTL,2); 

Thot = HTL(1,HTLsize);
Tcold = CTL(1,CTLsize);

temperatures = [CTL sampleTemps reversedHTL];



%Final calculations
TgradSample = (Thot-Tcold)/ySample; %Temperature gradient over the sample in K/m
kGraphite = AvgFlux/TgradSample; %Thermal conductivity of graphite sample in W/m*K
rGraphite = ySample/(kGraphite*areaSample); %Thermal resistance of graphite sample in K/W
rGraphiteAreaCompensated = rGraphite*areaSample; %Thermal resistance of graphite sample in K*m^2/W


%Analysis Code Ends

%Plotting Data:
figure(1);

%Plotting Temperatures vs Time
subplot(2,2,1);
plot(t,sH1);
hold on
plot(t,sH2);
hold on
plot(t,sH3);
hold on
plot(t,sC1);
hold on
plot(t,sC2);
hold on
plot(t,sC3);
hold on
xline(SSTime/1000,'--b',{'Steady', 'State'});
yticks(linspace(10,220,22));
ylim([10 220]);
grid on;
title(TestDataFilename,'Thermocouple Temperatures vs Time')
xlabel('Time (s)');
ylabel('Temperature (^{o}C)');


%Interpolation Plot
subplot(2,2,2);
%figure(2);
plot(length,HTL);
hold on
plot(length,HTL,'-x','MarkerIndices', 36000);
hold on
plot(length,HTL,'-s','MarkerIndices', [36000*0.122, 36000*0.5, 36000*0.877]);
hold on
plot(length,CTL);
hold on
plot(length,CTL,'-x','MarkerIndices', 36000);
hold on
plot(length,CTL,'-s','MarkerIndices', [36000*0.122, 36000*0.5, 36000*0.877]);
hold on
xline(0.036,'--b',{'Graphite', 'Surface', 'Temperature'})
yticks(linspace(10,220,22));
ylim([10 220]);
grid on;
title(TestDataFilename,'Interpolating Surface Temperature on Graphite Sample')
xlabel('Location on Al Block(m)');
ylabel('Temperature (^{o}C)');


%Temperatures Plot
subplot(2,2,3);
%figure(3);
plot(experimentalLength,temperatures);
hold on
plot(experimentalLength,temperatures,'-x','MarkerIndices', 36000);
hold on
plot(experimentalLength,temperatures,'-x','MarkerIndices', 36000 + ySample*10^6);
hold on
plot(experimentalLength, temperatures,'-s','MarkerIndices', [4392, 18000, 31572]);
hold on
plot(experimentalLength, temperatures,'-s','MarkerIndices', [36000 + ySample*10^6 + 4392, 36000 + ySample*10^6 + 18000, 36000 + ySample*10^6 + 31572]);
hold on
xline(0.036000,'--b'); %{'Graphite Cold Surface Temperature'}
hold on
xline((0.036000 + ySample),'--b'); %{'Graphite Hot Surface Temperature'}
yticks(linspace(10,220,22));
ylim([10 220]);
xlim([0 totalLength]);
grid on;
title(TestDataFilename,'Temperatures Across Test Setup');
xlabel('Location on Tester (m)');
ylabel('Temperature (^{o}C)');

pos = get(gcf,'Position');
set(gcf, 'Position',pos+[-500 -300 500 100])


%saveas(figure(1), fullfile('Graphs and Plots',TestDataFilename), 'png');


%Reporting Data
disp("For a Sample of thickness " + ySample + " m or " + ySample*10^3 + "mm" );
disp("Temperature Gradient over Hot Side is " + abs(HTslope) + " K/m");
disp("Temperature Gradient over Cold Side is " + CTslope + " K/m");
disp("Temperature Gradient over Sample is " + slopeSample + " K/m");
disp("Thermal Resistance of the Sample is " + rGraphiteAreaCompensated + " m^2 K/W");
%disp("Thermal Conductivity of the Sample is " + kGraphite + " W/ m K");




%---Function Definitions----------------------------------------------------




function valuesLinearFit = LinReg(Dataset)

    Distance = Dataset(:,2);
    Temperature = Dataset(:,1);
    
    %thickness_fit = (0:0.1:4)*1e-3;
    s = polyfit(Distance,Temperature,1);
    yfit = polyval(s,Distance);
    %fit_plot = polyval(s,thickness_fit);
    yres = Temperature-yfit;
    numerator = sum(yres.^2);
    denominator = (length(Temperature)-1)*var(Temperature);
    R2 = 1 - (numerator/denominator);
    
    m = s(1,1);
    c = s(1,2);
    
    
    valuesLinearFit = [ 
    m
    c
    R2
    ];

end

function plotValues = plotData(Fit,Dataset)

    slope = Fit(1,1);
    intercept = Fit(2,1);
    
    fitLine = slope.*Dataset(:,1) + intercept;
    
    plotValues = fitLine;

end
