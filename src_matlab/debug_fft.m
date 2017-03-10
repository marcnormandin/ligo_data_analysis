N = 10
data = zeros(N,1)
for n = 0:N-1
    data(n+1) = n - i*n;
end
data
y = ifft(data)
