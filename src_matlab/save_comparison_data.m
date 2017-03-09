% Save data to compare with the C version
formatSpec = '%f\t %f\t %f\n';
d = whitened_signal{1,1};
r = real(d);
im = imag(d);
fn = 'matlab_K1.whitened_signal';
disp(fprintf('Saving to whitened_signal to %s\n', fn))
disp(fn)
fid = fopen(fn, 'w');
for k = 1 : length(interp_fre)
    fprintf(fid, '%f\t %f\t %f\n', interp_fre(k), r(k), im(k));
end
fclose(fid);
