function [SFn]=HKI_Sub_SFn(SFnNum,eta,IntRange)
if IntRange(1) < 0 || IntRange(2) > 1
    SFn = 0;
elseif SFnNum == 1
    SFn = (eta-IntRange(1))/(IntRange(2)-IntRange(1)); % »ó½Â, x(num-1) = a < x < x(num) = b
elseif SFnNum == 2
    SFn = (IntRange(2)-eta)/(IntRange(2)-IntRange(1)); % ÇÏ°­, x(num) = b < x < x(num+1) = c
end
end