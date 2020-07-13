function  neuralData = problem_2(data, neuronNumber, plotVar)
 clc; close all;
out = pmdDataSetup(data);
reachAngles = out(:,1);
numTrials = length(reachAngles);
firingRates = out(:,2:end);
neuron = neuronNumber;
neuronFR = firingRates(:,neuron);
anglesC = cosd(reachAngles);
anglesS = sind(reachAngles);
 
offset = [anglesC,anglesS,ones(numTrials,1)];
coeff= regress(neuronFR,offset);
baseB = coeff(3);
k = sqrt(coeff(2)^2+coeff(1)^2);
pd = atand(coeff(2)/coeff(1));
if (plotVar)
    figure; 
    plot(reachAngles, neuronFR,'.');
    x = 0:350;
    y = baseB+k*cosd(pd-x);
    hold on 
    plot(x,y,'.');
    title('Cosine Wave Fitted to Tuning Curves')
    xlabel('Angle')
    ylabel('Firing Rates')
    hold off
end
neuralData = [baseB, k, pd];

