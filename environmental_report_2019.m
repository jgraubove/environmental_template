
clear all; clc; close all;

%--------------------------------------------------------%
% ISH Template for Environmental Data Analysis
% v.1.4
% Last modified: 14/10/2019
% Author: josep.grau.bove@ucl.ac.uk
%--------------------------------------------------------%
% Usage: This script reads the file "data.xlsx", which
% should be in the same folder.
%--------------------------------------------------------%
% Modifications & License:  This code is based on the 
% spreadsheet template created by N. Blades and extended 
% by K. Curran, M. Strlic and J. Grau over many years
% in the Institute for Sustainable Heritage. It has been 
% developed for teaching and research and is shared under 
% a GNU GPLv3 license (see "License" at the end of this file). 
% You are welcome to edit and improve this code and share 
% your improvements with the authors. 
%--------------------------------------------------------%


% Modifications on 15.10.2019 - n_days is added manually. 
%                             - the time axis is not automatic.

%Configuration
    
    save_plots = 1; % choose 0 to not save the plots 
                    % (it reduces computation time)
    
    % Environmental guidelines of choice (Default values for NT)
    disp([' '])
    minRH_ref = 40;
    maxRH_ref = 65;
    minT_ref = 5; 
    maxT_ref = 22;        

    % Read Data
    
    disp(['reading data...'])
    data = xlsread('data2.xlsx');


    % Allocate imported array to column variable names
    time = linspace(1,length(data(:,2)),length(data(:,2)));
    T = data(:,1);
    RH = data(:,2);
     
    % Number of data points in a day
    % This is used for the calculation of fluctuation limits
    
    n_day = 24; 

% (1) Humidity Calculations

    disp(['humidity calculations...'])


    for n = 1:length(T);
    
        % Vapour Saturation Pressure at reported T (in Pa)

        SatV(n) =10 ^ (28.59051 - (8.2 * (log10(T(n) + 273.16))) + (0.0024804 * (T(n) + 273.16)) -(3142.31 / (T(n) + 273.16)))*100000;

        % Vapour Saturation Pressuere at reported T and RH (in Pa)
        SatV2(n) =(RH(n)/100)*SatV(n);

        % Absolute Humidity (in g/kg)
        mabs(n) = 1000*0.622*(SatV2(n)/(101325-SatV2(n)));

        % Isoperm (Dimensionless)
        Rate_base = 50*exp(-100000/(8.314*(20+273)));
        I(n) = (RH(n)*exp(-100000/(8.314*(T(n)+273))))/Rate_base;
        
        
        
        % Lifetime Multiplier (Dimensionless)
        LM(n) = ((50/RH(n))^(1.3))*exp((-100000/8.314)*(1/(273+T(n))-1/293));

    end

% (2) Humidity Plots
    

    figure('Color',[1 1 1]); 
    yyaxis left
    plot(time/n_day,T); 
    ylabel('T(C)');
    
    yyaxis right
    plot(time/n_day,RH);
    ylabel('RH');
    xlabel('Time (days)'); ; 
    legend('T(C)','RH')
    
    if save_plots == 1; print('Figure 1 - TRH','-dpng', '-r300'); disp(['plots beig saved, it may take longer']); end;

    figure('Color',[1 1 1]); 
    plot(time/n_day,mabs,'b'); xlabel('Time (days)'); ; ylabel('Specific Humidity (g/kg)');title('Specific Humidity');
    if save_plots == 1; print('Figure 2 - Specific Humidity','-dpng', '-r300'); end;

    figure('Color',[1 1 1]); 
    plot(time/n_day,I,'g',time/n_day,LM); 
    xlabel('Time (days)'); ; 
    ylabel('Arbitrary units');
    title('Isoperm and Lifetime Multiplier');
    legend('Isoperm','Lifetime Multiplier');
    if save_plots == 1; print('Figure 3 - I','-dpng', '-r300'); end;

% (3) Comparison with Standards and Guidelines

    disp(['comparing with guidelines... '])

    figure('Color',[1 1 1]); 
    plot(RH,T,'o',-1,-1,'b',-1,-1,'g',-1,-1,'r',-1,-1,'m'); ylim([5 30]); xlim([10 80]);
    xlabel('RH'); ylabel('T(C)');
    grid MINOR;
        % Add Guidelines
        
        % NationalTrust
        
        minRH = 40;
        maxRH = 65;
        minT = 5; 
        maxT = 22; 
        
        r1 = rectangle('Position', [minRH,minT,maxRH-minRH,maxT-minT],'EdgeColor','b');
    
        % PD5454
        
        minRH = 35;
        maxRH = 60;
        minT = 13; 
        maxT = 20; 
        
        r1 = rectangle('Position', [minRH,minT,maxRH-minRH,maxT-minT],'EdgeColor','g');

        % ASHRAE AA / AICCM
        
        minRH = 45;
        maxRH = 55;
        minT = 15; 
        maxT = 25; 
        
        r1 = rectangle('Position', [minRH,minT,maxRH-minRH,maxT-minT],'EdgeColor','r');

        % V&A
        
        minRH = 40;
        maxRH = 65;
        minT = 18; 
        maxT = 25; 
        
        r1 = rectangle('Position', [minRH,minT,maxRH-minRH,maxT-minT],'EdgeColor','m');

    legend('Data','National Trust','PD5454','ASHRAE AA','V&A');
    if save_plots == 1; print('Figure 4 - TRH','-dpng', '-r300'); end;

    
    % Guideline Counter
    
    a = zeros(1,length(T)); b = a;
    
    for n=1:length(a); 
        if T(n)> maxT_ref;
            a(n) = a(n)+1;
        elseif T(n)< minT_ref;
            a(n) = a(n)+1;
        elseif RH(n)> maxRH_ref;
            b(n) = b(n)+1;
        elseif RH(n)< minRH_ref;
            b(n) = b(n)+1;
        end;
    end
    
    timeT = sum(a)/length(a);
    timeRH = sum(b)/length(b);
    

    
    
% (4) Mould Growth Calculation

    disp(['mould growth calculations...'])

    % Calculation of isopleth (Saildebauer single 
    Trm = linspace(0,30,100);

    for n=1:length(Trm) % Polynomial adjusted to data
    RHm(n) = -0.0017*Trm(n)^3 + 0.1144*Trm(n)^2 - 2.6784*Trm(n) + 97.017;
    end
    
    figure('Color',[1 1 1]); 
    plot(T,RH,'bl.',Trm,RHm,'r-','LineWidth',2); xlim([5 30]); ylim([10 100]);
    xlabel('T(C)'); ylabel('RH'); text(20,85,'Mould Germination');
    grid MINOR;
    if save_plots == 1; print('Figure 5 - Mould','-dpng', '-r300'); end;
    
    % Saildebauer regression coefficients
    
    COEFF = [-0.0012	-0.0012	-0.0012	-0.0012	-0.0012,
                    0.0914	0.0914	0.0914	0.0914	0.0914,
                    -2.37	-2.35	-2.36	-2.375	-2.36,
                    106.18	101	98	96	93];
    
    Ts = linspace(1,30,100); 
    
        for n=1:5; for i=1:100;
            Saildebauer(i,n) = COEFF(1,n)*Ts(i)^3+COEFF(2,n)*Ts(i)^2+COEFF(3,n)*Ts(i)+COEFF(4,n);
        end, end

    figure('Color',[1 1 1]); 
    plot(T,RH,'bl.',Ts,Saildebauer,'k-'); grid MINOR; ; xlim([0 30]); ylim([70 100]);
    xlabel('T(C)'); ylabel('RH'); 
        text(25,86,'0 days');
        text(25,81.5,'2 days');
        text(25,78,'4 days');
        text(25,75.5,'8 days');
        text(25,73,'16 days');


    
% (5) Fluctuation limits for hygroscopic materials (BS EN 15757)


    for n = 1:length(RH)-n_day*15-1;
        RH_average(n) = nan;
        sd(n) = nan;
    end
    for n = n_day*15:length(RH)-n_day*15-1;
        RH_average(n) = mean(RH(n+1-n_day*15:n+n_day*15));
        sd(n) = std(RH(n+1-n_day*15:n+n_day*15));
    end
    
    for n = length(RH)-n_day*15-1:length(RH);
        RH_average(n) = nan;
        sd(n) = nan;
    end
    
    figure('Color',[1 1 1]); 
    plot(time/n_day,RH,'Color',[0.7 0.7 0.7]); ylabel('RH'); xlabel('Time (days)'); ;;
    hold;
    plot(time/n_day,RH_average,'LineWidth',2)
    plot(time/n_day,RH_average+1.5*sd,'--r')
    plot(time/n_day,RH_average-1.5*sd,'--r')
    legend('Measured RH','Moving Average','Acceptable Fluctuations')
    if save_plots == 1; print('Figure 6 - Fluctuations','-dpng', '-r300'); end;

    
% (6) 2D histogram of T and RH data
    
    disp(['creating histogram...'])

    figure('Color',[1 1 1]); 
    histogram2(RH,T,50, 'YBinLimits', [5,30], 'XBinLimits', [10,80],...
    'DisplayStyle','tile',...
    'ShowEmptyBins','on',...
    'Normalization','probability');
    xlabel('RH'); ylabel('T (C)'); cb = colorbar; cb.Label.String = 'Normalised Probability';
    if save_plots == 1; print('Figure 6 - Histogram','-dpng', '-r300'); end;



% (7) Psychometric Chart

    disp(['creating pyschrometric chart...'])

    % NOTE. This chart is based on the "Simplified Psychrometric Chart"
    % by Izzi Urieli, downloadable at 
    % https://www.ohio.edu/mechanical/thermo/Applied/Chapt.7_11/Psychro_chart/psychro.html

    t = linspace(1,50,50); % temperature (C)
    % saturation vapor pressure (Pa)
    for n=1:length(t);
        pg(n) =10 ^ (28.59051 - (8.2 * (log10(t(n) + 273.16))) + (0.0024804 * (t(n) + 273.16)) -(3142.31 / (t(n) + 273.16)))*100000;
    end;
    patm = 101325; % standard atmosphere (Pa)
    rair = 0.287; % gas constant of air (kJ/kg.K)
    % saturation specific humidity
    for n=1:length(t);
        mg(n) = 1000*0.622*(pg(n)/(101325-pg(n)));
    end;

    figure('Color',[1 1 1]); 
    plot(t,mg,'r-')
    hold
    grid
    for phi = 0.1:0.1:0.4, % phi = relative humidity 10% - 40%
        m = 622*phi*pg./(patm-phi*pg);
        plot(t,m,'b')
    end
    for phi = 0.6:0.2:0.8, % phi = 60%, 80%
        m = 622*phi*pg./(patm-phi*pg);
        plot(t,m,'b')
    end

    % Text (surely there's a better way to include this)
    a=text(43,5,'10%');set(a,'Color','b');
    a=text(42,9,'20%');set(a,'Color','b');
    a=text(40,13,'30%');set(a,'Color','b');
    a=text(38,16,'40%');set(a,'Color','b');
    a=text(35,20,'60%');set(a,'Color','b');
    a=text(32,23,'80%');set(a,'Color','b');
    a=text(28,25,'100%');set(a,'Color','r');

    % specific volume and enthalpy/wet-bulb-temp
    t1 = [t(5),t(10),t(15),t(20),t(25),t(30),t(35)]; % saturation temperature (C)
    pg1 = [pg(5),pg(10),pg(15),pg(20),pg(25),pg(30),pg(35)]; % saturation pressure (kPa)
    wg1 = 622*pg1./(patm-pg1); % saturation specific humidity
    
    % specific volume of dry air (cubic m/kg dry air) 
    vol = rair.*(t1+273)./(patm-pg1); % specific vol at saturation
    tv0 = patm*vol/rair-273; % air temperature at zero humidity
    for i = 1:7,
        plot([t1(i),tv0(i)],[wg1(i),0],'g-')
    end

    a=text(28,27,'30 C');set(a,'Color','g');
    a=text(24,21,'25 C');set(a,'Color','g');
    a=text(18,15,'20 C');set(a,'Color','g');
    a=text(14,11,'15 C');set(a,'Color','g');
    a=text(8,8,'10 C');set(a,'Color','g');
    a=text(4,6,'5 C');set(a,'Color','g');


    % Uncomment the following to plot Enthalpy 

%         % wet bulb temperature (also enthalpy) lines (red)
%         h = t1 + 2.5*wg1 % enthalpy (kJ/kg-dry-air) (displayed)
%         t0 = h; % temperature at zero humidity for enthalpy h
%         for i = 1:6,
%         	plot([t1(i),t0(i)],[wg1(i),0],'r-')
%         end
%         % enthalpy axis and enthalpy lines (black)
%         for h = 10:10:110, % enthalpy (kJ/kg-dry-air)
%         	t0 = h; % temperature at zero humidity
%         	t1 = (h - 12.5)/3.5; % temperature on the enthalpy axis
%         	w1 = t1 + 5; % specific humidity on the enthalpy axis
%         	plot([t0,t1],[0,w1],'k-')
%         end
%         plot([0,25],[5,30],'k-') % the oblique enthalpy axis
        
        
    grid minor;
    scatter(T,mabs,'.');
    axis([0,50,0,30]) % limit the range of the chart
    title('Psychrometric Chart')
    xlabel('Dry Bulb Temperature (C)')
    ylabel('Specific Humidity (g/kg)')
        if save_plots == 1; print('Figure 7 - Chart','-dpng', '-r300'); end;

% (8) Equilibrium Moisture Content (Experimental)
    
    disp(['emc calculations... '])

    % Calculation of T averaged 24 hours 
    
    
    for n = 1:length(T)-n_day*0.5-1;
        T_average(n) = nan;
        sd(n) = nan;
    end
    for n = n_day*15:length(T)-n_day*0.5-1;
        T_average(n) = mean(T(n+1-n_day*0.5:n+n_day*0.5));
        sd(n) = std(T(n+1-n_day*0.5:n+n_day*0.5));
    end
    
    for n = length(T)-n_day*0.5-1:length(T);
        T_average(n) = nan;
        sd(n) = nan;
    end
    

    
    for n=1:length(T)
        
        RHi = RH_average/100;...
        Ti = (T_average*(9/5))+32;
                  
        W(n) = 330 + 0.452*Ti(n) + 0.00415*Ti(n)^2;
        K(n) = 0.791 + 0.000463*Ti(n) - 0.000000844*Ti(n)^2;
        K1(n) = 6.34 + 0.000775*Ti(n) - 0.0000935*Ti(n)^2;
        K2(n) = 1.09 + 0.0284*Ti(n) - 0.0000904*Ti(n)^2;
        M(n) = (1800/W(n))*((K(n)*RHi(n)/(1-K(n)*RHi(n))) + ((K1(n)*K(n)*RHi(n) + (2*K1(n)*K2(n)*K(n)^2*RHi(n)^2)) / (1 + K1(n)*K(n)*RHi(n) + K1(n)*K2(n)*K(n)^2*RHi(n)^2)));

    end
    
    figure('Color',[1 1 1]); 
    plot(time/n_day,M,'r',[max(time/n_day),min(time/n_day)],[7,7],...
        [max(time/n_day),min(time/n_day)],[10.5,10.5],...
        [max(time/n_day),min(time/n_day)],[12.5,12.5],'k--',...
        [max(time/n_day),min(time/n_day)],[5,5],'k--');
    legend('%EMC','Acceptable Risk for Metals (IPI)','High Risk for Metals (IPI)','Safe region for Hygroscopic Materials')


    xlabel('Time (days)'); ; ylabel('EMC (%)'); title('Equilbrium Moisture Content');
    ylim([0,20])
    if save_plots == 1; print('Figure 8 - EMC','-dpng', '-r300'); end;
    
% (Final) Screen displays
    disp([' '])
    disp(['RESULTS'])
    disp([' '])
    disp(['Min.RH ', num2str(min(RH)),' Max.RH ', num2str(max(RH)), ' Range ', num2str(max(RH)-min(RH))])
    disp(['Min.T ', num2str(min(T)),' Max.T ', num2str(max(T)), ' Range ', num2str(max(T)-min(T))])
    disp([' '])
    disp(['T is in specification for ', num2str((1-timeT)*100),' % of the time'])
    disp(['RH is in specification for ', num2str((1-timeRH)*100),' % of the time'])
    disp([' '])

%--------------------------------------------------------%
%  License:
%--------------------------------------------------------%
% Copyright (C) 2017  Josep Grau-Bove
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% This code can be run with zero or minamal changes 
% in Octave, Scilab and Matlab. The user is however
% advised that the use of Matlab-specific instructions may render 
% the code unsuitable for the GNU license, as the related libraries
% are not necessarily free software. 
% 
% --------------------------------------------------------%
   

