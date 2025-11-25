fs = 44100
t = 1:(fs*10);
f = 1;
A = 0*t;
while f<(fs/2)
  f = f * exp(1/137.036);
  A = A + sin(2*pi*f*t/fs);
end
soundsc(A,fs,16);
