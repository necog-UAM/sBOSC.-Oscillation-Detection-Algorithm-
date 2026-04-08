function [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2)

if isempty(c0s) %%%% ly correction
    A = NaN;
    mu = NaN;
    sigma = NaN;
    rsq = NaN;

else

    parfor i = 1:length(c0s)
        c0 = c0s(i);

        x = centers;
        y = counts.*rectp(find(foi2==c0),:);    % isolate individual peaks

        f = fit(x.' ,y.','gauss1');

        A(i)  = f.a1;
        mu(i) = f.b1;
        sigma(i) = f.c1/sqrt(2);

        gausf = A(i) * exp(-(x-mu(i)).^2 /(2*sigma(i).^2));

        yresid = counts - gausf;
        SSresid = sum(yresid.^2);
        SStotal = (length(counts)-1) * var(counts);
        rsq(i) = 1 - (SSresid/SStotal) ;
    end
end

