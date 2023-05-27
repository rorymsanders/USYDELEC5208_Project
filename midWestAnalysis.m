%the purpose of this script is to analyse the optimisation results. Data is
%pulled in from the results of the optimisation and categorised into the
%various loss classes
clear;
clc;
%pull data from excel file
SolarMX = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"C2:C745");
WindMX = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"D2:D745");
LoadMX = readmatrix('ProjectDataTimeseries.xlsx','Sheet','InputDataHour','Range',"E2:N745");

%ResultsMX = readmatrix('ProjectDataTimeseries.xlsx','Sheet','Results','Range',"B2:ABQ34");
ResultsMX = readmatrix('ProjectDataTimeseries.xlsx','Sheet','ResultsNoBess','Range',"B2:ABQ34");

ResultsMX = transpose(ResultsMX);

%add a total column to the LoadMX
for i = 1:17520
    LoadMX(i,11) = sum(LoadMX(i,1:10));
end

%initialise and indicator matrix for the results
IndiMx = zeros(17520,17);

%define the limits of the transmissions lines
Limits = ...
    [71.3
    32.4;
    34.8;
    30.1;
    21;
    77.5;
    27;
    16.92;
    32.6;
    27;
    27;
    27;];

%pull in the ResultsMX and flag if one of the 12 lines has hit a limit by
%comparing the lines load flow to it's limit. Store the loadflow headroom
%in the IndiMx
for i = 1:17520
    
   for j = 1:12
       if abs(ResultsMX(i,j+3))==Limits(j)
           IndiMx(i,4) = 1;
           IndiMx(i,5:16) = transpose(Limits) - abs(ResultsMX(i,4:15));
       end
   end
    
    
end

i = 0;

for i = 1:17520
    %check if possible power is larger than actual generation and no line limits. Indicates constrained
    %energy due to surplus renewable power
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 0
        IndiMx(i,1) = 1;
    else
        IndiMx(i,1) = 0;
    end
    
    %check if actual from generator 1 and 2 is less than possible power and
    %at least one line limit has been reached but no grid supply
    %Indicates constrained energy due to line limits
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 1 && ResultsMX(i,3)==0
        IndiMx(i,2) = 1;
    else
        IndiMx(i,2) = 0;    
    end
    
    %check if actual from generator 1 and 2 is less than possible power.
    %at least one line limit requiring grid power purchasing
    %Indicates constrained energy due to line limits
    if (SolarMX(i,1)+WindMX(i,1))>sum(ResultsMX(i,1:2)) && IndiMx(i,4) == 1 && ResultsMX(i,3)>0
        IndiMx(i,3) = 1;
    else
        IndiMx(i,3) = 0;    
    end
    
end

i = 0;

sum(IndiMx(:,1));
sum(IndiMx(:,2));
sum(IndiMx(:,3));
sum(IndiMx(:,4));

totalLoad = sum(LoadMX(:,11)*0.5);
totalPossible = sum(SolarMX(:,1)+WindMX(:,1))*0.5;
totalRenewableGen = sum(sum(ResultsMX(:,1:2)*0.5));
totalGridGen = sum(ResultsMX(:,3)*0.5);

totalWastedExcess = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(:,1)+ResultsMX(:,2)))*0.5.*IndiMx(:,1));
totalWastedLineLimitExcess = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(:,1)+ResultsMX(:,2)))*0.5.*IndiMx(:,2));
totalWastedLineLimitCurtailment = sum(((SolarMX(:,1)+WindMX(:,1))-(ResultsMX(:,1)+ResultsMX(:,2)))*0.5.*IndiMx(:,3));

totalPossible-totalRenewableGen
totalWastedExcess+totalWastedLineLimitExcess+totalWastedLineLimitCurtailment

        