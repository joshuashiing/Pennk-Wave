function fr = recons_conj(f)

% Based on the conjugate relation, reconstruct the spectrum
N = length(f);

if rem(N, 2) == 0
    nh = (N - 2) / 2;
    fp = f(2 : (nh+1));
    fn = flip(conj(fp));
    fr = [real(f(1)); fp(:); real(f(nh+2)); fn(:)];
else
    nh = (N - 1) / 2;
    fp = f(2 : (nh+1));
    fn = flip(conj(fp));
    fr = [real(f(1)); fp(:); fn(:)];
end
fr = reshape(fr, size(f));