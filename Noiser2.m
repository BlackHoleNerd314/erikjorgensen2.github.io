midi = -24;
fs = 44100;
f = 440*2^((midi-9)/12);
t = 0:(1/fs):ceil(440/f);
A = (440/f);
A0 = exp(-f*t/24);
y = A0.*erf(A*sin(2*pi*f*t));
sound(y,fs,16);




