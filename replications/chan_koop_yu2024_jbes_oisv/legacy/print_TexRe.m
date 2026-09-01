% This function transforms the outputs into LaTex version outputs (with stars).

function re = print_Texre(v, star, option)
re = [];
if option == 1 [~, idx] = min(v); elseif option == 2 [~, idx] = max(v); end
for i = 1:length(v)
    vi = num2str(v(i),'%.3f');
    stari = star(i); word = '';
    if stari == 1 word = '\ensuremath{^{*}}';
    elseif stari == 2 word = '\ensuremath{^{**}}';
    elseif stari == 3 word = '\ensuremath{^{***}}';
    end
    if i == idx
        re_str = strcat('\textbf{',vi,'}',word,' ');
    else
        re_str = strcat(vi,word,' ');
    end
    re = [re; split(re_str)];
end
