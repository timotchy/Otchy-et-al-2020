for i = 1:length(p)
    temp2 = sum(p{i}-q{i});
    stretchy(i) = temp2;
end
barh(reverse(stretchy))
ylim([1 length(p)])