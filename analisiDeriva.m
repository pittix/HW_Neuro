function [media,varianza] = analisiDeriva(ROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
numSamples =  size(ROI,1);
slope = zeros(numSamples,1);
intercept = zeros(numSamples,1);
time = 1:1:size(ROI(1).tac,1);
for i=1:1: numSamples
  [slope(i),intercept(i)]= fit(  time, ROI(1).tac,'poly1' );
  fitParabola = fit(  time, ROI(1).tac,'poly2' );
end

media = mean(slope);
varianza = var(slope);
end

