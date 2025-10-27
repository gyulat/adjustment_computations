% draw samples from the Tukey supermodel
function x = draw_tukey(p,sc,nsamples)
    component = rand(1, nsamples) < (1-p);
    x = zeros(1, nsamples);
    x(component)  = 0 + 1 * randn(1, sum(component));
    x(~component) = 0 + sc * randn(1, sum(~component));
end
