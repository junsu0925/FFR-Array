function IntResult = HKI_Sub_Integral(Int_SFn,IntRange,IntError)
Method = 'auto';
IntResult=integral2(Int_SFn,IntRange(1)+IntError,IntRange(2)-IntError,0,2*pi,'Method',Method,'AbsTol',1e-12,'RelTol',1e-8);
end