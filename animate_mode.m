
%mode=201; %gravest gravity/Kelvin mode = 201.
%mode=204;
if ~exist('mode', 'var')
    mode = numstrm + 1;
end

figure(2); clf; colormap(darkjet);

% Get angular frequency,
freq = myeigs(mode);

% Get period, so we'll plot one period.
T = 2*pi/abs(freq);
numframes = 50;
t=0:T/numframes:T;  %make array of times to evaluate mode at.

for jj=1:length(t)
    real_eta = real(exp(1i*freq*t(jj))*eta{mode});
    if jj == 1
        eta_max = max(abs(real_eta(:)));
    end
    pf2d(N, x, y, real_eta); colorbar;
    view([0 90]); 
    axis tight; axis equal;
    caxis([-eta_max eta_max]);
    colormap(darkjet);
    title(['w/f=' sprintf('%1.2f, t=%1.2f of %1.2f',scaledeigs(mode),t(jj),T)]);
    print('-dpng',sprintf('frame%07d.png',jj));
    drawnow;
end
