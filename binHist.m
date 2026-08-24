binsWanted = [5 6];

function h = binHist(m, binsWanted, name, col)
    d  = hypot(m.tx - m.apex(1), m.ty - m.apex(2));
    in = false(size(d));
    for b = binsWanted
        in = in | (d >= m.edges(b) & d < m.edges(b+1));
    end
    h = m.coh(in);
    fprintf('%-8s bins %s: n %4d  mean %.3f  p90 %.3f  frac>0.6 %.3f\n', ...
        name, mat2str(binsWanted), numel(h), mean(h), prctile(h,90), mean(h>0.6));
end

hA = binHist(m_aci, binsWanted, 'ACI', 1);
hB = binHist(m_bn,  binsWanted, 'BN',  2);

figure('Color','w');
histogram(hA, 0:0.05:1, 'Normalization','probability', ...
          'FaceColor',[0.85 0.35 0.2], 'FaceAlpha',0.55); hold on
histogram(hB, 0:0.05:1, 'Normalization','probability', ...
          'FaceColor',[0.2 0.4 0.75], 'FaceAlpha',0.55);
legend('ACI','BN'); xlabel('local coherence'); ylabel('fraction of windows');
title(sprintf('bins %s', mat2str(binsWanted)));

[~, pKS] = kstest2(hA, hB);
fprintf('\nKS test on the two distributions: p = %.3g\n', pKS);