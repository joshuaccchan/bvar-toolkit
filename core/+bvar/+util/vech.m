% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace). Canonical source: chan2023_joe_mlvarsv/legacy/utility/vech.m (single legacy copy).
function y=vech(Y)
y = nonzeros(tril(Y));
end
