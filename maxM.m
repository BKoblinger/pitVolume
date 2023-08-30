function [ maxVal, r, c ] = maxM( A )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[maxVal, location] = max(A(:));

[r,c] = ind2sub(size(A),location); 

end

