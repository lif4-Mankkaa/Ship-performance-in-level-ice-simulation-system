function [ grouped, noOfContacts ] = grouping( vectorin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vectorin = vectorin(:);

vectorin_m = vertcat(-1, vectorin, vectorin(end) + 2);
indexDiff = diff(vectorin_m);
len = diff(find(indexDiff ~= 1));
noOfContacts = length(len);
grouped = zeros(1,2 * noOfContacts);
start = 1;
j = 1;

for i = 1:2:2 * noOfContacts
    grouped(i) = vectorin(start);
    if len(j) ~= 1
        grouped(i+1) = vectorin(len(j)-1+start);
        start = start + len(j);
    else
        grouped(i+1) = vectorin(start);
        start = start + 1;
    end
    j = j + 1;
end

end

