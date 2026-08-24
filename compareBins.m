function comparesBins(m_aci, m_bn, binsWanted)

   binsWanted = [5 6];

dA = hypot(m_aci.tx - m_aci.apex(1), m_aci.ty - m_aci.apex(2));
inA = false(size(dA));
for b = binsWanted, inA = inA | (dA >= m_aci.edges(b) & dA < m_aci.edges(b+1)); end
hA = m_aci.coh(inA);

dB = hypot(m_bn.tx - m_bn.apex(1), m_bn.ty - m_bn.apex(2));
inB = false(size(dB));
for b = binsWanted, inB = inB | (dB >= m_bn.edges(b) & dB < m_bn.edges(b+1)); end
hB = m_bn.coh(inB);

fprintf('ACI n %4d  mean %.3f  p90 %.3f  frac>0.6 %.3f\n', ...
        numel(hA), mean(hA), prctile(hA,90), mean(hA>0.6));
fprintf('BN  n %4d  mean %.3f  p90 %.3f  frac>0.6 %.3f\n', ...
        numel(hB), mean(hB), prctile(hB,90), mean(hB>0.6));

figure('Color','w');
histogram(hA, 0:0.05:1, 'Normalization','probability', ...
          'FaceColor',[0.85 0.35 0.2], 'FaceAlpha',0.55); hold on
histogram(hB, 0:0.05:1, 'Normalization','probability', ...
          'FaceColor',[0.2 0.4 0.75], 'FaceAlpha',0.55);
legend('ACI','BN'); xlabel('local coherence'); ylabel('fraction of windows');