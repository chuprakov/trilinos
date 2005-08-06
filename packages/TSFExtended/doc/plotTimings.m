
A = load('vectorTimings.dat')

n = A(:,1);

tOp = A(:,3);

tCore = A(:,5);

plot(log10(n), tCore./tOp, '-o') ; axis([0 6 0 1.1]) ; line([0 6], [1 1])

xlabel('log (N)')

ylabel('TSF overloaded operator timings / Thyra timings')

title('efficiency of operator overloading relative to core computations')