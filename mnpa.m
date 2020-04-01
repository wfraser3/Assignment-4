R1 = 1;
cap = 0.25;
R2 = 2;
L = 0.2;
R3 = 7.8091; %From findR3Faster.m
alpha = 100;
R4 = 0.1;
R0 = 1000;

%{ 
Unknowns:

V1
V2
V3
V4
V0
IL

%}

G = [1 0 0 0 0 0;
    1/R1 ((-1/R1)-(1/R2)) 0 0 0 -1;
    0 0 -1/R3 0 0 1;
    0 0 -alpha/R3 1 0 0;
    0 0 (-alpha*R0)/R3 0 (1+R0) 0;
    0 -1 1 0 0 0];

C = [0 0 0 0 0 0;
    cap -cap 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 L];

%a) DC Case
counter = 1;
vout = zeros(21,1);
v3 = zeros(21,1);
for loop = -10:10
    F = [loop 0 0 0 0 0];
    V = G\F';
    vout(counter) = V(5);
    v3(counter) = V(3);
    counter = counter + 1;
end

vin = linspace(-10,10,21);

figure(1)
plot(vin,vout,'-*b')
xlabel('Vin (V)')
ylabel('Vout (V)')
title('Relationship between Vin and Vout (DC Case)')
xlim([-10,10])
ylim([min(vout),max(vout)])

figure(2)
plot(vin,v3,'-*b')
xlabel('Vin (V)')
ylabel('V3 (V)')
title('Relationship between Vin and V3 (DC Case)')
xlim([-10,10])
ylim([min(v3),max(v3)])

%b) AC Case
w = linspace(0,1000,1001);
f = w/(2*pi);
voutAC = zeros(1001,1);
F = [1 0 0 0 0 0];
for loop = 1:1001
acMatrix = G + (1i*w(loop)*C);
V = acMatrix\F';
voutAC(loop) = V(5); 
end

magV = abs(voutAC);

figure(3)
plot(f,magV,'b')
xlabel('Frequency (Hz)')
ylabel('Vout (V)')
title('AC Case')

gain = 20*log10(magV);

figure(4)
plot(f,gain,'b')
xlabel('Frequency (Hz)')
ylabel('Gain (dB)')
title('AC Gain')

%c) Pertubations of C
miu = cap;
std = 0.05;
w = pi;
cVector = linspace(0.01,0.5,50);
cProb = (1/(std*sqrt(2*pi)))*exp(-0.5*(((cVector-miu)./std).^2));
cProb = cProb*10;

start = 1;
probVector = zeros(round(sum(cProb)),1);
for loop = 1:length(cProb)
    prob = round(cProb(loop));
    if(prob==1)
        probVector(start) = cVector(loop);
        start = start + 1;
    elseif(prob>1)
        stop = start + prob - 1;
        probVector(start:stop) = cVector(loop);
        start = stop + 1;
    end
end    

probVector = probVector(randperm(length(probVector)));
clear gain V voutAC magV
gain = zeros(1000,1);

for time = 1:1000
    index = randperm(length(probVector),1);
    newCap = probVector(index);
    
    C = [0 0 0 0 0 0;
    newCap -newCap 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 L];
    
    acMatrix = G + (1i*w*C);
    V = acMatrix\F';
    voutAC = V(5); 
    magV = abs(voutAC);
    gain(time) = 20*log10(magV);
    
    
    figure(5)
    histogram(gain(1:time),50,'FaceColor','b')
    xlabel('Gain (dB)')
    ylabel('Counts')
    title('Histogram of Gain with Pertubations of C')
    
    pause(0.01)
    
end

%Transient
figureCounter = 6;
Gtran = [1 0 0 0 0 0;
    1/R1 ((-1/R1)-(1/R2)) 0 0 0 -1;
    0 0 -1/R3 0 0 1;
    0 0 -alpha/R3 1 0 0;
    0 0 (-alpha*R0)/R3 0 (1+R0) 0;
    0 -1 1 0 0 0];

Ctran = [0 0 0 0 0 0;
    cap -cap 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 L];

acMatrix = Ctran + Gtran;

time = zeros(1001,1);
for i = 2:1001
    time(i) = (i-1)/1000;
end

TDtitle{1} = 'Time Domain Transient Response to Step Function';
TDtitle{2} = 'Time Domain Transient Response to Sine Function';
TDtitle{3} = 'Time Domain Transient Response to Gaussian Function';

FDtitle{1} = 'Frequency Domain Transient Response to Step Function';
FDtitle{2} = 'Frequency Domain Transient Response to Sine Function';
FDtitle{3} = 'Frequency Domain Transient Response to Gaussian Function';

livePlotting = input('Enter 1 for live plotting, 0 for no live plotting: ');

for part = 1:3
    if(part == 1)
        vin = zeros(1001,1);
        vin(61:1001) = 1;
    elseif(part == 2)
        f = 1/0.03;
        vin = sin(2*pi*f*time);
    elseif(part == 3)
        b = 0.06;
        c = 0.03;
        vin = exp((-(time - b).^2)/(2*(c^2)));        
    end

    for loop = 1:1001
        F = [vin(loop) 0 0 0 0 0];
        currentTime = (loop-1);
        if(loop==1)
            Vtran = acMatrix\F';
            lastVtran = Vtran;
            lastVoutTran = abs(Vtran(5));
            saveVout{part}(loop) = abs(Vtran(5));
        elseif(loop==2)
            tempF = F + (Ctran*lastVtran);
            Vtran = acMatrix\tempF';
            saveVout{part}(loop) = abs(Vtran(5));
            voutTran = abs(Vtran(5));
            if(livePlotting==1)
                figure(figureCounter)
                plot([time(loop-1),time(loop)],[vin(loop-1),vin(loop)],'b','DisplayName','Vin')
                hold on
                plot([time(loop-1),time(loop)],[lastVoutTran,voutTran],'g','DisplayName','Vout')
                xlabel('Time (s)')
                ylabel('Voltage (V)')
                title(['Time = ',num2str(currentTime),'ms'])
                figureCounter = figureCounter + 1;
            end
            lastVoutTran = voutTran;
            lastVtran = Vtran;
        else
            tempF = F + (Ctran*lastVtran);
            Vtran = acMatrix\tempF';
            saveVout{part}(loop) = abs(Vtran(5));
            voutTran = abs(Vtran(5));
            if(livePlotting==1)
                figure(figureCounter-1)
                plot([time(loop-1),time(loop)],[vin(loop-1),vin(loop)],'b','DisplayName','Vin')
                plot([time(loop-1),time(loop)],[lastVoutTran,voutTran],'g','DisplayName','Vout')
                title(['Time = ',num2str(round(currentTime)),'ms'])
            end
            lastVoutTran = voutTran;
            lastVtran = Vtran;
        end
        if(livePlotting==1)
            pause(0.001)
        end
    end
    
    if(livePlotting==0)
        figure(figureCounter)
        plot(time,vin,'b','DisplayName','Vin')
        hold on
        plot(time,saveVout{part},'g','DisplayName','Vout')
        title(TDtitle{part})
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        legend
        hold off
        figureCounter = figureCounter + 1;
    end
    
    fs = 1000;
    n = 1000;
    f = (-n/2:n/2-1)*(fs/n);
    fourierIn{part} = fft(vin);
    fourierIn{part} = fftshift(fourierIn{part});
    fourierIn{part} = abs(fourierIn{part});
    fourierOut{part} = fft(saveVout{part});
    fourierOut{part} = fftshift(fourierOut{part});
    fourierOut{part} = abs(fourierOut{part});
    
    figure(figureCounter)
    plot(f,fourierIn{part}(1:end-1),'b','DisplayName','Input')
    hold on
    plot(f,fourierOut{part}(1:end-1),'g','DisplayName','Output')
    set(gca,'FontSize',18)
    title(FDtitle{part})
    xlabel('Frequency (Hz)')
    ylabel('|FFT|')
    legend
    figureCounter = figureCounter + 1;
end