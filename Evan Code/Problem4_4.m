
%

[spikes, pos, vel] = Problem4_0_1(Data);

% I would index to the size of the spikingData variable. INDEXING BY THE
% SIZE OF SPIKES SMALLEST MATRIX

% firingRateMean=mean(spikes,1);
% firingRateCovariance=cov(spikes);
% 
% II = firingRateMean; Mx1
% 
% V = firingRateCovariance; MxM
% stateModel =

% Zt = A*Zt(i) + Noise # i refers to the previous trial
% input Zt is a velocity value that will be used to get an approximate
% postion value

% STEPS
% 1) Preprocessing -> initiating the A, C, Q, R variables with appropriate
% size MxM, use pmdDataSetup, make sure that there are 50 of each 
% 2) Initializing first time stamp, 

%% SORTING DATA BY REACH ANGLE
data = pmdDataSetup(Data);

angles = [0 45 90 135 180 225 270 315];
reachAngles = data(:,1);
count = 0;
trainingData.veloX = {808};
trainingData.veloY = {808};
trainingData.posX = {808};
trainingData.posY = {808};
trainingData.spikes= {808};

sortedIndex = [];
for i =1:length(angles)
    index = reachAngles==angles(i);
    for j=1:length(index)
        if index(j)==1&&count<101
            sortedIndex = [sortedIndex j];
            count = count +1;
        elseif index(j)==1&&count>=101
            count = 0;
            break;
        end
    end
end


%%
for i = 1:length(sortedIndex)
    trainingData.veloX{i} = vel.x{sortedIndex(i)};
    vel.x{sortedIndex(i)}(:) = [];
    trainingData.veloY{i} = vel.y{sortedIndex(i)};
    vel.y{sortedIndex(i)}(:) = [];
    trainingData.posX{i} = pos.x{sortedIndex(i)};
    pos.x{sortedIndex(i)}(:) = [];
    trainingData.posY{i} = pos.y{sortedIndex(i)};
    pos.y{sortedIndex(i)}(:) = [];
    trainingData.spikes{i} = spikes{sortedIndex(i)};
    spikes{sortedIndex(i)} = [];
end



  
  
  %% A & C UPDATE
  
  firstASum = zeros(4,4);
  secondASum = zeros(4,4);
  firstSumC = zeros(231,4);
  secondSumC = zeros(4,1);
  Xt = zeros(231,1);

  
   for trial=1:size(trainingData.spikes) % outer loop goes through all the trails 
       
       T=size(trainingData.spikes{trial}, 1); % inner loop iterates through bin sizes of spike counts per trial
        
        % redefine Zt for each time bin 

        for binPos=2:T 
            
         Xt = trainingData.spikes{trial}(binPos,:); % D = numNeurons needs to be 233 x 1

         Xt = Xt';
            
          if binPos == T + 1
       
              continue;
          end 
          

        if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
        end 
        
         if(isnan(trainingData.veloX{trial}(binPos - 1)) || isnan(trainingData.veloY{trial}(binPos - 1)) || isnan(trainingData.posX{trial}(binPos - 1)) || isnan(trainingData.posY{trial}(binPos-1)))
            
            continue;
            
        end 
        
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
        
        % check isnan for vel/pos values, if so then continue 
           
       velX2 = trainingData.veloX{trial}(binPos - 1);
       velY2 = trainingData.veloY{trial}(binPos - 1);
       posX2 = trainingData.posX{trial}(binPos - 1);
       posY2 = trainingData.posY{trial}(binPos - 1);
        
        
        Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
        
        Zt2 = [posX2, posY2, velX2, velY2]';
                
        firstASum =  (Zt*Zt2') + firstASum;

        secondASum = (Zt2*Zt2') + secondASum;
        
        firstSumC = (Xt*Zt') + firstSumC;
        
        secondSumC = (Zt*Zt') + secondSumC;
        
        end                
        
   end 
   %% COMPUTE A & C

     
   A = firstASum * (secondASum)^(-1); % is MxM so 4x4, computing from only Z1 and Z2
   C = firstSumC * (secondSumC)^(-1); % D x M multiplying a 233x1 by a 1x4 results in 233 x 4

 
 
   %% Q & R UPDATE
   
   QSum = zeros(4,4);
   firstTermQ = 0;
   RSum = zeros(231,231);
   Xt = zeros(231,1);
   firstTermR = 0;
   
    for trial=1:size(trainingData.spikes) % outer loop goes through all the trails 
         
       T=size(trainingData.spikes{trial}, 1); % inner loop iterates through bin sizes of spike counts per trial
        
        % redefine Zt for each time bin 

        for binPos=2:T 
            
             Xt = trainingData.spikes{trial}(binPos,:)'; % D = numNeurons needs to be 233 x 1
             
             if binPos == T + 1
       
              continue;
             end 
             
         if(isnan(trainingData.veloX{trial}(binPos)) || isnan(trainingData.veloY{trial}(binPos)) || isnan(trainingData.posX{trial}(binPos)) || isnan(trainingData.posY{trial}(binPos)))
            
            continue;
            
        end 
        
         if(isnan(trainingData.veloX{trial}(binPos - 1)) || isnan(trainingData.veloY{trial}(binPos - 1)) || isnan(trainingData.posX{trial}(binPos - 1)) || isnan(trainingData.posY{trial}(binPos-1)))
            
            continue;
            
        end 
          
       velX1 = trainingData.veloX{trial}(binPos);
       velY1 = trainingData.veloY{trial}(binPos);
       posX1 = trainingData.posX{trial}(binPos);
       posY1 = trainingData.posY{trial}(binPos);
        
       velX2 = trainingData.veloX{trial}(binPos - 1);
       velY2 = trainingData.veloY{trial}(binPos - 1);
       posX2 = trainingData.posX{trial}(binPos - 1);
       posY2 = trainingData.posY{trial}(binPos - 1);
        
        
        
        Zt = [posX1, posY1, velX1, velY1]'; %changed order to pos, vel
        
        Zt2 = [posX2, posY2, velX2, velY2]';
            
        
        QSum = (Zt-A*Zt2)*((Zt-A*Zt2)') + QSum;
       
        RSum = (Xt-C*Zt)*((Xt-C*Zt)') + RSum;
        
        end
        
        firstTermQ = (T-1)+ firstTermQ;
        
        firstTermR = (T + firstTermR);
    end 
      
   %% COMPUTE Q & R
   
  Q = (1/firstTermQ) * QSum; % M x M
  R = (1/firstTermR) * RSum;
%   R = (1/firstTermR) * RSum; % D x D (233 x 233)

%% Testing Phase
%% PREPROCESSING THE TESTING DATA, removing empty cells, put in new structs

for i=1:size(pos.x,2)
    if isempty(pos.x{i})
        pos.x{i} = 0;
    end 
     if isempty(pos.y{i})
        pos.y{i} = 0;
     end  
     if isempty(vel.x{i})
        vel.x{i} = 0;
     end 
     if isempty(vel.y{i})
        vel.y{i} = 0;
     end 
     if isempty(spikes{i})
        spikes{i} = 0;
    end 
end

posXRemovedZeros = pos.x(cellfun(@(x) ~isequal(x, 0), pos.x));
posYRemovedZeros = pos.y(cellfun(@(x) ~isequal(x, 0), pos.y));
velXRemovedZeros = vel.x(cellfun(@(x) ~isequal(x, 0), vel.x));
velYRemovedZeros = vel.y(cellfun(@(x) ~isequal(x, 0), vel.y));
spikesRemovedZeros = spikes(cellfun(@(x) ~isequal(x, 0), spikes));



%% LOOPING KALMAN FUNCTION

predictedValues.muPos(length(posXRemovedZeros)) = {1};
predictedValues.muVelo(length(posXRemovedZeros)) ={1};
predictedValues.sig(length(posXRemovedZeros)) = {1};

Z_1_values = [];


for trial=1:size(spikesRemovedZeros,2) % outer loop goes through all the trails 
         
       T=size(spikesRemovedZeros{trial}, 1); % inner loop iterates through bin sizes of spike counts per trial
        
       velX1 = velXRemovedZeros{trial};
       velY1 = velYRemovedZeros{trial};
       posX1 = posXRemovedZeros{trial};
       posY1 = posYRemovedZeros{trial};
    
       Z_1_values = [[posX1(1), posY1(1), velX1(1), velY1(1)]; Z_1_values];
       
       Z_1_mean = mean(Z_1_values)';
       Z_1_covar = cov(Z_1_values);
end 
%%

for trial=1:size(spikesRemovedZeros,2) % outer loop goes through all the trails 

     T=size(spikesRemovedZeros{trial}, 1); % inner loop iterates through bin sizes of spike counts per trial

       for binPos = 2:T
           
       mu = zeros(size(Z_1_mean)); % 2 by bin size should be 4x1
       sig = zeros(size(Z_1_covar,2));
       mu1 = zeros(size(Z_1_mean,2))'; % 4x1
       sig1 = zeros(size(Z_1_covar,2)); % 4x4
           
       Xt = spikesRemovedZeros{trial}(binPos,:); % D = numNeurons needs to be 231 x 1

       Xt = Xt';
       
       if binPos == 1
           
        mu1 = Z_1_mean;
        sig1 = Z_1_covar;
        
        else 
            
        mu1 = A*mu;
        sig1 = A*sig*A' + Q;
        
       end 
       
        Kt = sig1*C'*(C*sig1*C'+R); % dimension of 231x4
        mu = mu1 + Kt*(Xt - C*mu1); % dimension of 4x1
        sig = sig1 - Kt*C*sig1; 
        
        if binPos == 1
             predictedValues.muPos{trial} = mu(1:2)';
             predictedValues.muVelo{trial} =mu(3:4)';
             predictedValues.sig{trial} = sig;
            
        else 
        
           predictedValues.muPos{trial} = [predictedValues.muPos{trial}; mu(1:2)'];
           predictedValues.muVelo{trial} = [predictedValues.muVelo{trial}; mu(3:4)'];       
           predictedValues.sig{trial} = [predictedValues.sig{trial}; sig];
       
        end 
       end
end 


figure;
for trial =1:length(posXRemovedZeros)
    if trial==1
         plot(predictedValues.muPos{trial}(:,1), predictedValues.muPos{trial}(:,2), 'Color', 'b')
        plot(posXRemovedZeros{trial}, posYRemovedZeros{trial})
        hold on
    else
         plot(predictedValues.muPos{trial}(:,1), predictedValues.muPos{trial}(:,2), 'Color', 'k')
         plot(posXRemovedZeros{trial}, posYRemovedZeros{trial})
    end
end









