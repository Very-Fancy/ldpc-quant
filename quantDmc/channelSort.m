function [Pout,ind] = channelSort(Pin)
%channelSort - Sort the outputs of a binary-input DMC.
%
%   With DMC PIN, where rows sum to 1, POUT = CHANNELSORT(PIN)
%   gives the sorted channel.
%
%   PIN is a DMC conditional distribution, and POUT is sorted to 
%   satisfy the condition:
%
%      Pr(Y=1|X=1)     Pr(Y=2|X=1)             Pr(Y=M|X=1)
%     ------------- < ------------- <  ...  < -------------
%      Pr(Y=1|X=2)     Pr(Y=2|X=2)             Pr(Y=M|X=2)
%
%    which corresponds to
%
%        POUT(1,1) / POUT(2,1)  < POUT(1,2) / POUT(2,2)  < ... < POUT(1,M) / POUT(2,M)
%
%   for a 2-by-M matrix PIN, POUT.
%
%   Additionally, POUT is normalized so each row sums to one.  A random 
%   channel can be generated by:
%
%     CHANNELSORT(RAND(2,M))
%
% QuantDMC (c) Brian Kurkoski and contributors
% Distributed under an MIT-like license; see the file LICENSE

[J M] = size(Pin);
assert(J==2,'Pin must have two rows')

Pin(1,:) = Pin(1,:) / sum(Pin(1,:));
Pin(2,:) = Pin(2,:) / sum(Pin(2,:));

[~,ind] = sort( Pin(1,:) ./ Pin(2,:) );
Pout = Pin(:,ind);
