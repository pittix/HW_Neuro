function [media,varianza] = analisiDeriva(ROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
numSamples =  size(ROI,1);
slope = zeros(numSamples,1);
intercept = zeros(numSamples,1);
time = (1:1:size(ROI,2))';

for i=1:1: numSamples
  fitObj= fit(time,  ROI(i,:)','poly1' );
  slope(i) = fitObj.p1;
  intercept(i) = fitObj.p2;
end

media = mean(slope);
varianza = var(slope);
end

