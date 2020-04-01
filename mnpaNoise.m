R1 = 1;
cap = 0.25;
R2 = 2;
L = 0.2;
R3 = 7.8091; %From findR3Faster.m
alpha = 100;
R4 = 0.1;
R0 = 1000;
Cn = 0.00001;

%{ 
Unknowns:

V1
V2
V3
V4
V0
IL

%}

figureCounter = 1;

%Transient
miu = 0.1;
std = 0.01;
w = pi;
iVector = linspace(miu-(3*std),miu+(3*std),1000);
iProb = (1/(std*sqrt(2*pi)))*exp(-0.5*(((iVector-miu)./std).^2));
iProb = iProb*10;

start = 1;
probVector = zeros(round(sum(iProb)),1);
for loop = 1:length(iProb)
    prob = round(iProb(loop));
    if(prob==1)
        probVector(start) = iVector(loop);
        start = start + 1;
    elseif(prob>1)
        stop = start + prob - 1;
        probVector(start:stop) = iVector(loop);
        start = stop + 1;
    end
end    

probVector = probVector(randperm(length(probVector)));

%figureCounter = 6;
Gtran = [1 0 0 0 0 0;
    1/R1 ((-1/R1)-(1/R2)) 0 0 0 -1;
    0 0 -1/R3 0 0 1;
    0 0 -alpha/R3 1 0 0;
    0 0 (-alpha*R0)/R3 0 (1+R0) 0;
    0 -1 1 0 0 0];

Ctran = [0 0 0 0 0 0;
    cap -cap 0 0 0 0;
    0 0 -Cn 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 L];

acMatrix = Ctran + Gtran;

time = zeros(1001,1);
for i = 2:1001
    time(i) = (i-1)/1000;
end

TDtitle{3} = 'Time Domain Transient Response to Gaussian Function with Noise';

FDtitle{3} = 'Frequency Domain Transient Response to Gaussian Function with Noise';

livePlotting = input('Enter 1 for live plotting, 0 for no live plotting: ');

part = 3;

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
    index = randperm(length(probVector),1);
    In(loop) = probVector(index);
    %In(loop) = rand;
    F = [vin(loop) 0 In(loop) 0 0 0];
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