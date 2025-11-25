N = 10
a0 = zeros(N+1,N+1);
for n=1:N
  a0(n,n+1)=sqrt(n);
end
a0;
a1 = [0,1;0,0];
a1;
a0*a0-a0*a0;
a0'*a0'-a0'*a0';
a0*a0'-a0'*a0;
a1*a1+a1*a1;
a1'*a1'+a1'*a1';
a1*a1'+a1'*a1;

a = kron(eye(N),a1);
A = kron(eye(N),a1');

ans = a*A+A*a







