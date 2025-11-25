function [data] = testDerv0(m0,r)
  m = exp(randn(1)*1)*m0;
  a = exp(randn(1)*1)*m0;
  x = exp(randn(1)*1)*r;
  y = exp(randn(1)*1)*r;
  z = exp(randn(1)*1)*r;
  err = sqrt(sum(sum(abs(KerrMetric(m,a,x,y,z)-Metric(m,a,x,y,z)).^2)));
  err0 = sqrt(sum(sum(sum(abs((Connection(m,a,x,y,z)-Connection0(m,a,x,y,z))).^2))));
  data = [err,err0]
end
