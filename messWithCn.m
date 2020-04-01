R1 = 1;
cap = 0.25;
R2 = 2;
L = 0.2;
R3 = 7.8091; %From findR3Faster.m
alpha = 100;
R4 = 0.1;
R0 = 1000;
Cn = [0.00001,1,10,100];

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

for cloop = 1:length(Cn)
    %Transient
    miu = 0.5;
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
    
    Gtran = [1 0 0 0 0 0;
        1/R1 ((-1/R1)-(1/R2)) 0 0 0 -1;
        0 0 -1/R3 0 0 1;
        0 0 -alpha/R3 1 0 0;
        0 0 (-alpha*R0)/R3 0 (1+R0) 0;
        0 -1 1 0 0 0];

    Ctran = [0 0 0 0 0 0;
        cap -cap 0 0 0 0;
        0 0 -Cn(cloop) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 L];

    acMatrix = Ctran + Gtran;

    time = zeros(1001,1);
    for i = 2:1001
        time(i) = (i-1)/1000;
    end

    FDtitle{3} = 'Noise Bandwidth with Different Values of Cn';

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
            lastVoutTran = voutTran;
            lastVtran = Vtran;
        else
            tempF = F + (Ctran*lastVtran);
            Vtran = acMatrix\tempF';
            saveVout{part}(loop) = abs(Vtran(5));
            voutTran = abs(Vtran(5));
            lastVoutTran = voutTran;
            lastVtran = Vtran;
        end
    end

    fs = 1000;
    n = 1000;
    f = (-n/2:n/2-1)*(fs/n);
    fourierIn{cloop} = fft(vin);
    fourierIn{cloop} = fftshift(fourierIn{cloop});
    fourierIn{cloop} = abs(fourierIn{cloop});
    fourierOut{cloop} = fft(saveVout{part});
    fourierOut{cloop} = fftshift(fourierOut{cloop});
    fourierOut{cloop} = abs(fourierOut{cloop});

    figure(figureCounter)
    hold on
    plot(f(451:551),fourierOut{cloop}(451:551),'DisplayName',['C = ',num2str(Cn(cloop))])
    set(gca,'FontSize',18)
    title(FDtitle{part})
    xlabel('Frequency (Hz)')
    ylabel('|FFT|')
    legend
end