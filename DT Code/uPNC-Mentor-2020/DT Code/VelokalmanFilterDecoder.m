function predictedValues = VelokalmanFilterDecoder(data, veloData, spikes)

D = length(spikes{1}(2,:));
M = 2;

%%
testingCoeff = 0.90;
testingsize = floor(length(data)*testingCoeff);
while rem(testingsize,8)~=0
    testingsize = testingsize +1;
end
%%
angles = [0 45 90 135 180 225 270 315];
reachAngles = data(:,1);
count = 0;
trainingData.veloX = {testingsize};
trainingData.veloY = {testingsize};
trainingData.spikes= {testingsize};

sortedIndex = [];
for i =1:length(angles)
    index = reachAngles==angles(i);
    for f=1:length(index)
        if index(f)==1&&count<(testingsize/8)
            sortedIndex = [sortedIndex f];
            count = count +1;
        elseif index(f)==1&&count>=(testingsize/8)
            count = 0;
            break;
        end
    end
end
%%
for i = 1:length(sortedIndex)
    trainingData.veloX{i} = veloData.x{sortedIndex(i)};
    veloData.x{sortedIndex(i)}= [];
    veloData.x{sortedIndex(i)} = 0;
    trainingData.veloY{i} = veloData.y{sortedIndex(i)};
    veloData.y{sortedIndex(i)} = [];
    veloData.y{sortedIndex(i)} = 0;
    trainingData.spikes{i} = spikes{sortedIndex(i)};
    spikes{sortedIndex(i)} = [];
    spikes{sortedIndex(i)} = 0;
end
%%
testingVeloX = veloData.x(cellfun(@(j) ~isequal(j,0), veloData.x));
testingVeloY = veloData.y(cellfun(@(j) ~isequal(j,0), veloData.y));
testingSpikes = spikes(cellfun(@(j) ~isequal(j,0), spikes));


%%

QSum= zeros(M,M);
QSum1= zeros(M,M);
firstTermQ= 0;
firstTermR= 0;
Rsum = zeros(D,D);
%Prediction 
for trial =1:length(trainingData.spikes)
    
    trialveloX = trainingData.veloX{trial};
    trialveloY = trainingData.veloY{trial};
    trialSpikes = trainingData.spikes{trial};  
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end
                if bin~=1&&(isnan(trialveloX(bin-1))||isnan(trialveloY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 
                end
                if bin ==1
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [veloX2 veloY2]';
                    Ctop = (x_t'*z_t');
                    Cbottom = (z_t*z_t');
                elseif bin ==2
                    veloX = trialveloX(bin-1);
                    veloY = trialveloY(bin-1);
                    z_t1 = [veloX veloY]';
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [veloX2 veloY2]';
                    Atop =(z_t*z_t1');
                    Abottom =(z_t1*z_t1');
                    Ctop = Ctop+ (x_t'*z_t');
                    Cbottom = Cbottom+ (z_t*z_t');
                else
                    veloX = trialveloX(bin-1);
                    veloY = trialveloY(bin-1);
                    z_t1 = [veloX veloY]';
                    veloX2 = trialveloX(bin);
                    veloY2 = trialveloY(bin);
                    x_t = trialSpikes(bin,:);
                    z_t = [veloX2 veloY2]';
                    Atop = Atop+ (z_t*z_t1');
                    Abottom = Abottom+ (z_t1*z_t1');
                    Ctop = Ctop +(x_t'*z_t');
                    Cbottom = Cbottom+ (z_t*z_t');
                end
        end
end
A = Atop*Abottom^(-1);
C = Ctop*Cbottom^(-1);
for trial =1:length(trainingData.spikes)
    trialveloX = trainingData.veloX{trial};
    trialveloY = trainingData.veloY{trial};
    trialSpikes = trainingData.spikes{trial};
    T = size(trialSpikes,1);
        for bin =1:T
                if (isnan(trialveloX(bin))||isnan(trialveloY(bin))||isnan(trialSpikes(bin)))
                    continue;
                end

                if bin~=1&&(isnan(trialveloX(bin-1))||isnan(trialveloY(bin-1))||isnan(trialSpikes(bin-1)))
                    continue; 

                end
                    if bin==1
                        veloX = trialveloX(bin);
                        veloY = trialveloY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [veloX veloY]';
                        R = (x_t' - C*z_t)*(x_t'-C*z_t)';
                        Rsum = R + Rsum;
   
                    else
                        veloX = trialveloX(bin-1);
                        veloY = trialveloY(bin-1);
                        z_t1 = [veloX veloY]';
                        veloX2 = trialveloX(bin);
                        veloY2 = trialveloY(bin);
                        x_t = trialSpikes(bin,:);
                        z_t = [veloX2 veloY2]';
                        QSum1 = (z_t-A*z_t1)*(z_t-A*z_t1)' +QSum1; %+ QSum1;
                        R =(x_t' - C*z_t)*(x_t'-C*z_t)';
                        Rsum = R + Rsum;
                    end
                    firstTermQ = (T-1)+ firstTermQ;
                    firstTermR = (T)+ firstTermQ;
        end
end

Q= (1/firstTermQ) * QSum1;
R = (1/firstTermR) *Rsum;
%% testing phase
predictedValues.muVelo(length(testingVeloX)) ={1};
predictedValues.sig(length(testingVeloX)) = {1};
Z_1 = zeros(length(testingVeloX),M);

for trial=1:length(testingVeloX)
    Z_1(trial,:) = [testingVeloX{trial}(1) testingVeloY{trial}(1)];
end
Z_1_mean = mean(Z_1)';
Z_1_covar = cov(Z_1);

%%
for trial=1:length(testingVeloX)
    T =size(testingSpikes{trial},1);
    Mu_t1 = zeros(size(Z_1_mean));
    Mu_t = Mu_t1;
    Sig_t = zeros(size(Z_1_covar));
    Sig_t1 = Sig_t;
    for bin =1:T
       if bin ==1
           Mu_t1 = Z_1_mean;
           Sig_t1 = Z_1_covar;

       else
           Mu_t1 = A*Mu_t;
           Sig_t1 = A*Sig_t*A'+Q;
       end
       k_t = Sig_t1*C'*(C*Sig_t1*C'+R)^(-1);
       Mu_t= Mu_t1+k_t*(testingSpikes{trial}(bin,:)'-C*Mu_t1);
       Sig_t = Sig_t1-k_t*C*Sig_t1;
       if bin==1
           predictedValues.muVelo{trial} =Mu_t(1:2)';
           predictedValues.sig{trial} = Sig_t;
       else 
           predictedValues.muVelo{trial} = [predictedValues.muVelo{trial}; Mu_t(1:2)'];
           predictedValues.sig{trial} = [predictedValues.sig{trial}; Sig_t];     
       end
    end 
end
%% PREFORMANCE CALCULATION

% Postion & Velocity Difference

error.diffXVelo(length(testingVeloX)) = {1};
error.diffYVelo(length(testingVeloX)) = {1};
meanErrorVeloX = zeros(length(testingVeloX),1);
meanErrorVeloY = zeros(length(testingVeloX),1);

for trial= 1:length(testingVeloX)
    
       if (max(max(isnan(predictedValues.muVelo{trial}(:,1)))) == 1 || max(max(isnan(predictedValues.muVelo{trial}(:,2)))) == 1 || max(max(isnan(testingVeloX{trial}(1:end-1)))) == 1 || max(max(isnan(testingVeloY{trial}(1:end-1)))) == 1)
           continue;
       end
       xv = abs(predictedValues.muVelo{trial}(:,1) - testingVeloX{trial}(1:end-1))/(testingVeloX{trial}(1:end-1));
       yv = abs(predictedValues.muVelo{trial}(:,2) - testingVeloY{trial}(1:end-1))/(testingVeloY{trial}(1:end-1));
       error.diffXVelo{trial} = xv((xv ~= 0));
       error.diffYVelo{trial} = yv((yv ~= 0));
       meanErrorVeloX(trial) = mean(error.diffXVelo{trial},1);
       meanErrorVeloY(trial) = mean(error.diffYVelo{trial},1);

end
predictedValues.Errorperformance= mean([mean(abs(meanErrorVeloX))*100, mean(abs(meanErrorVeloY))*100]);
