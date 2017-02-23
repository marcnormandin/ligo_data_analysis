cd /Users/marcnormandin/build/bin

fid = fopen('K1')
C = textscan(fid, '%f')

plot(C{1,1})