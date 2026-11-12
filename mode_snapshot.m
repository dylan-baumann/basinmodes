function s = mode_snapshot(e)
% MODE_SNAPSHOT  Real-valued snapshot of a complex mode field, at peak amplitude.
%
%   s = mode_snapshot(eta{k})
%
%   basinmodes.m returns each mode as a COMPLEX field (eta{k}, potential{k},
%   streamfcn{k}), because the underlying eigenproblem is complex. The overall
%   complex phase of an eigenvector is arbitrary: eig() fixes it by a
%   normalisation convention that has nothing to do with the physics. So
%   real(eta{k}) is "the mode at a random instant in its cycle", and for a
%   standing seiche that instant can be the moment the water is flat
%   (real(eta) ~ 0), giving a blank or noise-scaled picture.
%   (For a nearly non-rotating basin this happens for almost every mode;
%   with rotation it happens for some of them.)
%
%   This function multiplies by exp(-i*phi) so that the point of largest
%   |eta| is real and positive, i.e. it returns the snapshot at the instant of
%   maximum displacement. For a standing wave the result is the real mode
%   shape; for a travelling/rotating wave it is the instant when the largest
%   displacement is at its crest.
%
%   Input and output have the same size (Np x K).

    [~, imax] = max(abs(e(:)));
    s = real(e * exp(-1i * angle(e(imax))));
end
