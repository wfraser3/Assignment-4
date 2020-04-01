%{
ELEC 4700: Assignment 2
William Fraser
101001393
%}

%{
        L - (i) 
     W   ___________________________
     |  |1 8                        |
    (j) |2 ?                        |
        |3                          |
        |4                          |
        |5                          | 
        |6                          | 
        |7__________________________|
%}

W = 20;
L = 40;

cmap = ones(W,L);
contactsCond = 10^(-2);
cmap(1:8,15:25) = contactsCond;
cmap(12:20,15:25) = contactsCond;

G = zeros(W*L);
B = zeros(1,W*L);

voltSweep = linspace(0.1,10,10);
xCurrentDensity = zeros(10,1);

for sweeper = 1:10
    display(['Voltage = ',num2str(voltSweep(sweeper)),'V'])
    for i = 1:L
        for j = 1:W
            n = j + (i-1)*W;
            if(i==1)
                G(n,:) = 0;
                G(n,n) = 1;
                B(n) = voltSweep(sweeper);
            elseif(i==L)
                G(n,:) = 0;
                G(n,n) = 1;
            elseif(j==1)
                nxm = j + (i-2)*W;
                nxp = j + i*W;
                nyp = j + 1 + (i-1)*W;

                resXM = (cmap(j,i)+cmap(j,i-1))/2;
                resXP = (cmap(j,i)+cmap(j,i+1))/2;
                resYP = (cmap(j,i)+cmap(j+1,i))/2;

                G(n,n) = -(resXM+resXP+resYP);
                G(n,nxm) = resXM;
                G(n,nxp) = resXP;
                G(n,nyp) = resYP;
            elseif(j==W)
                nxm = j + (i-2)*W;
                nxp = j + i*W;
                nym = j-1 + (i-1)*W;

                resXM = (cmap(j,i)+cmap(j,i-1))/2;
                resXP = (cmap(j,i)+cmap(j,i+1))/2;
                resYM = (cmap(j,i)+cmap(j-1,i))/2;

                G(n,n) = -(resXM+resXP+resYM);
                G(n,nxm) = resXM;
                G(n,nxp) = resXP;
                G(n,nym) = resYM; 
            else
                nxm = j + (i-2)*W;
                nxp = j + i*W;
                nym = j-1 + (i-1)*W;
                nyp = j + 1 + (i-1)*W;

                resXM = (cmap(j,i)+cmap(j,i-1))/2;
                resXP = (cmap(j,i)+cmap(j,i+1))/2;
                resYM = (cmap(j,i)+cmap(j-1,i))/2;
                resYP = (cmap(j,i)+cmap(j+1,i))/2;

                G(n,n) = -(resXM+resXP+resYM+resYP);
                G(n,nxm) = resXM;
                G(n,nxp) = resXP;
                G(n,nyp) = resYM;
                G(n,nym) = resYP; 
            end
        end
    end

    T = G\B';

    voltage = zeros(W,L);
    row = 1;
    column = 1;

    for loop = 1:length(T)
        voltage(row,column) = T(loop);
        if(row == W)
            row = 1;
            column = column + 1;
        else
            row = row + 1;
        end
    end

    for j = 1:L
        for i = 1:W
            if i == 1
                dropY(i, j) = (voltage(i + 1, j) - voltage(i, j))/5;
            elseif i == W
                dropY(i, j) = (voltage(i, j) - voltage(i - 1, j))/5;
            else
                dropY(i, j) = (voltage(i + 1, j) - voltage(i - 1, j))/10;
            end
            if j == 1
                dropX(i, j) = (voltage(i, j + 1) - voltage(i, j))/5;
            elseif j == L
                dropX(i, j) = (voltage(i, j) - voltage(i, j - 1))/5;
            else
                dropX(i, j) = (voltage(i, j + 1) - voltage(i, j - 1))/10;
            end
        end
    end

    dropX = -dropX;
    dropY = -dropY;

    x = zeros(L*W,1);
    y = zeros(L*W,1);
    index = 1;

    for i = 1:W
        for j = 1:L
            x(index) = j;
            y(index) = i;
            Evector(index,1) = dropX(i,j);
            Evector(index,2) = dropY(i,j);
            Cvector(index,1) = cmap(i,j);
            index = index + 1;
        end
    end

    Xcurrent = Cvector.*Evector(:,1);
    Ycurrent = Cvector.*Evector(:,2);
    
    xCurrentDensity(sweeper) = mean(Xcurrent);
end

current = xCurrentDensity*100;
figure(1)
plot(current,voltSweep,'-*b')
xlabel('Current (A)')
ylabel('Voltage (V)')
title('Current - Voltage Characteristic of the Device')

R3 = (voltSweep(10)-voltSweep(1))/(current(10)-current(1));