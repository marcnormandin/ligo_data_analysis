cd /Users/marcnormandin/build/bin

detId = cellstr(['H1';'L1';'V1';'K1']);
f = figure
for id = 1:1:length(detId)
    compare_plot(id, detId, interp_fre, whitened_signal)
end

function compare_plot(id, detId, interp_fre, whitened_signal)
subplot(2,2,id)
fid = fopen(detId{id})
C = textscan(fid, '%f %f')
fclose(fid)
plot(C{1,1}, C{1,2}, 'r')
xlabel('Frequency')
ylabel('Real part of chirp')
%title('Comparison between Matlab and C versions')
title(detId{id})
hold all
plot(interp_fre, real(whitened_signal{1,id}), 'b:')
legend('C version', 'Matlab version')
end

