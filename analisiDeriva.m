function [slope] = analisiDeriva(ROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
TR=2.6;
numSamples =  size(ROI,1);
slope = zeros(numSamples,1);
time = (1:1:size(ROI,2))';

for i=1:1: numSamples
  fitObj= fit(TR*time,  ROI(i,:)','poly1' );
  slope(i) = fitObj.p1;
end


end

