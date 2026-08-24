function varargout = radialOrientation(varargin)
% radialOrientation  Vessel alignment relative to the optic nerve head.
%
%   s = radialOrientation(skel, apexXY, bandMask, pxUM, minLenUM)
%   radialOrientation('selftest')      % <-- run this first
%
% Returns, for the branches whose centroid lies inside bandMask:
%   radialAlignment   +1 all radial, -1 all circumferential, 0 no preference
%   radialCoherence    0 relative angle scattered, 1 relative angle consistent
%   globalCoherence    0 no common absolute direction, 1 all parallel
%   meanRadialAngle_deg   0 = outward, 90 = circumferential, 45 = random
%   nUsed, meanLen_um
%
% ---------------------------------------------------------------------
% MATH
%   Branch axis v: unit eigenvector of the larger eigenvalue of the
%     covariance of the branch's (x,y) pixel coordinates. Computed here
%     rather than taken from regionprops 'Orientation', so that v and the
%     radial vector are guaranteed to live in the same coordinate frame.
%   Radial unit vector r = (centroid - apex) / |centroid - apex|.
%   cos(delta) = |v . r|   - absolute because a vessel is an undirected
%     line, so delta is folded into [0, 90] degrees.
%   radialAlignment = sum(w * cos(2 delta)) / sum(w),  cos2d = 2cos^2 - 1
%   radialCoherence = |sum(w * exp(2i delta))| / sum(w)
%   globalCoherence = |sum(w * exp(2i theta))| / sum(w),
%     theta = atan2(vy, vx). Angle doubling is the standard treatment for
%     axial data: it makes 0 and 180 degrees identical.
%   w = branch length in pixels, so long segments count more than short.
% ---------------------------------------------------------------------

    if nargin == 1 && ischar(varargin{1}) && strcmpi(varargin{1}, 'selftest')
        selftest();
        return
    end
    [varargout{1:nargout}] = measure(varargin{:});
end

% ======================================================================
function s = measure(skel, apexXY, bandMask, pxUM, minLenUM)

    if nargin < 5 || isempty(minLenUM), minLenUM = 20; end
    minPix = max(round(minLenUM / pxUM), 5);

    nbr    = conv2(double(skel), [1 1 1;1 0 1;1 1 1], 'same');
    juncPx = skel & nbr >= 3;
    L      = bwlabel(skel & ~juncPx, 8);

    st = regionprops(L, 'PixelList', 'Centroid');   % PixelList is [x y]

    cosd_  = []; theta = []; w = [];

    for k = 1:numel(st)
        P = st(k).PixelList;
        if size(P,1) < minPix, continue; end

        c = st(k).Centroid;                          % [x y]
        if bandMask(round(c(2)), round(c(1))) == 0, continue; end

        C = cov(double(P));                          % 2x2 over [x y]
        [V, Dg] = eig(C);
        [~, im] = max(diag(Dg));
        v = V(:, im);  v = v / norm(v);              % branch axis, [vx; vy]

        r = [c(1) - apexXY(1); c(2) - apexXY(2)];
        nr = norm(r);
        if nr < eps, continue; end
        r = r / nr;                                  % radial unit vector

        cosd_(end+1,1) = abs(v' * r);                %#ok<AGROW>  cos(delta)
        theta(end+1,1) = atan2(v(2), v(1));          %#ok<AGROW>
        w(end+1,1)     = size(P,1);                  %#ok<AGROW>
    end

    if isempty(w)
        s = struct('radialAlignment',NaN,'radialCoherence',NaN, ...
                   'globalCoherence',NaN,'meanRadialAngle_deg',NaN, ...
                   'nUsed',0,'meanLen_um',NaN);
        return
    end

    W     = sum(w);
    delta = acos(min(max(cosd_,0),1));               % [0, pi/2]

    s.radialAlignment    = sum(w .* (2*cosd_.^2 - 1)) / W;
    s.radialCoherence    = abs(sum(w .* exp(2i*delta))) / W;
    s.globalCoherence    = abs(sum(w .* exp(2i*theta))) / W;
    s.meanRadialAngle_deg = sum(w .* rad2deg(delta)) / W;
    s.nUsed              = numel(w);
    s.meanLen_um         = mean(w) * pxUM;
end

% ======================================================================
function selftest()
% Four patterns with known answers. If these do not come out as expected,
% do not trust the measurement on real data.

    N = 601; c = (N+1)/2; apex = [c c];
    band = true(N);
    pxUM = 1;

    fprintf('\n%-22s %10s %10s %10s %10s %6s\n', 'pattern', ...
        'radAlign', 'radCoher', 'globCoher', 'meanAng', 'n');
    fprintf('%s\n', repmat('-', 1, 74));

    % --- 1. spokes radiating from the centre -> radial ---
    BW = false(N);
    for a = 0:15:359
        for rr = 40:1:280
            x = round(c + rr*cosd(a)); y = round(c + rr*sind(a));
            BW(y, x) = true;
        end
    end
    show('spokes (expect +1, 0)', measure(BW, apex, band, pxUM, 20));

    % --- 2. concentric rings -> circumferential ---
    BW = false(N);
    for rr = 60:40:280
        for a = 0:0.2:359
            x = round(c + rr*cosd(a)); y = round(c + rr*sind(a));
            BW(y, x) = true;
        end
    end
    show('rings (expect -1, 0)', measure(BW, apex, band, pxUM, 20));

    % --- 3. parallel horizontal lines -> globally aligned, not radial ---
        BW = false(N);
    for y = 40:40:N-40
        for x0 = 60:80:N-140
            BW(y, x0:x0+60) = true;
        end
    end
    show('parallel (expect 0, 1)', measure(BW, apex, band, pxUM, 20));

    % --- 4. randomly oriented segments -> no preference either way ---
    rng(1);
    BW = false(N);
    for k = 1:400
        a  = rand*180;  len = 30;
        x0 = 60 + rand*(N-120); y0 = 60 + rand*(N-120);
        for t = 0:0.5:len
            x = round(x0 + t*cosd(a)); y = round(y0 + t*sind(a));
            if x>=1 && x<=N && y>=1 && y<=N, BW(y,x) = true; end
        end
    end
    show('random (expect 0, 0)', measure(BW, apex, band, pxUM, 20));

    fprintf(['\nradAlign  +1 radial, -1 circumferential, 0 none\n' ...
             'globCoher  1 all parallel to one axis, 0 none\n' ...
             'Spokes and rings are the key pair: both are radially\n' ...
             'symmetric so globCoher cannot separate them, but\n' ...
             'radAlign should be near +1 and -1 respectively.\n\n']);
end

function show(name, s)
    fprintf('%-22s %10.3f %10.3f %10.3f %10.1f %6d\n', name, ...
        s.radialAlignment, s.radialCoherence, s.globalCoherence, ...
        s.meanRadialAngle_deg, s.nUsed);
end