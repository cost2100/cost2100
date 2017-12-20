function test_poissrnd_own
% TEST_POISSRND_OWN. As its name indicates, this routine tests the
% poissonness of the  random number generator `poissrnd_own'. Values of
% lambda in {1,2,5,10,20,50,100} are tested. 
% 
% For each such value 200 samples with 500 elements each are considered.
% The chi-square goodness-of-fit test metric, V, is computed for each
% sample. The KS test is then applied to the resulting 200 values of V, and
% the samples are declared as Po(lambda) distributed (H0), or not (H1), at
% 5 percent significance. 
%
% For the sake of comparison, the test is also applied to `poissrnd', in
% Matlab's Statistics Toolbox.
%
% As a EXAMPLE, type
%
% > clear all; close all; clc; profile clear; profile on; tic; test_poissrnd_own; toc; profile viewer;
%
% and check the result of the KS test (in the legend of the generated
% plots.)
%
% The CONCLUSION is that based on the above test, we cannot reject the null
% hypothesis. Indeed, for the values of lambda examined, `poissrnd_own' is
% as good as `poissrnd'. Matlab's implementation runs faster; but that
% seems due to vectorization. A vectorized version `poissrnd_own_parallel'
% is also provided (but not used, since `poissrnd_own' is not called very
% often). The following table provides a comparison of the runtime of the
% two algorithms.
%
%         Avg. serial runtime [s]   Avg. runtime (parallel) [s]
% Matlab          126.0                          2.6
% Own              16.7                          2.2
%
% The original algorithm was published by D.E. Knuth in [1]. A more
% accurate version for large values of lambda appears in [2]. See also [3]
% for additional Poisson number generators. For the general theory of
% hypothesis testing, see [4]. Detailed information about the chi-squared
% goodness-of-fit test can be found in [5]--[7]. Finally, for a concise
% overview of both the Chi-squared and the Kolmogorov-Smirnov test, the
% reader is referred to [1].
% 
% References:
% [1] `The art of computer programming, vol. 2 (3rd ed.): seminumerical
% algorithms', D.E. Knuth.
% [2] https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
% [3] https://en.wikipedia.org/wiki/Poisson_distribution#Confidence_interval
% [4] https://onlinecourses.science.psu.edu/statprogram/node/136
% [5] http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm
% [6] http://www.stat.yale.edu/Courses/1997-98/101/chigf.htm
% [7] `The Chi-squared test of goodness of fit', William G. Cochran, Ann.
% of Math. Stat., Vol 23, No. 3 (Sep. 1952), pp. 315--345.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This file is a part of the COST2100 channel model.
%
%This program, the COST2100 channel model, is free software: you can 
%redistribute it and/or modify it under the terms of the GNU General Public 
%License as published by the Free Software Foundation, either version 3 of 
%the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but 
%WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
%or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
%for more details.
%
%If you use it for scientific purposes, please consider citing it in line 
%with the description in the Readme-file, where you also can find the 
%contributors.
%
%You should have received a copy of the GNU General Public License along 
%with this program. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set to one if you want to test the vectorized implementations (Matlab's
% versus own)
VECTORIZATION_ON = 1;

% Sample size.
n = 500;

% Sample size and binning. We want to run the chi-square goodness-of-fit
% test with at least 5 counts in each bin. Using a sample size n=500, we
% obtain the following usable bins:
lambda_list = [ 1   2   5   10    20    30    50   100 ].';
bin_lo_list = [ 0   0   1    4    13    22    41    90 ].';
bin_hi_list = [ 3   4   9   15    26    37    58   110 ].';

for lambda_idx = 1:length(lambda_list)
    
    % For each targetted intensity, run a large number of tests   
    lambda = lambda_list(lambda_idx);
    bin_lo = bin_lo_list(lambda_idx);
    bin_hi = bin_hi_list(lambda_idx);
    nbins = bin_hi - bin_lo + 1;
    
    fprintf('lambda = %d\n',lambda);

    % In the chi-square test, the number of degrees of freedrom is the number
    % of bins minus 1
    dof = nbins - 1;

    V_ref = zeros(200,1);
    V_own = zeros(200,1);

    for k = 1:200
        
        % Obtain the samples and classify them into bins
        expCounts = zeros(nbins,1);
        expCounts(1) = sum(n * (lambda.^(0:bin_lo).' * exp(-lambda) ./ factorial((0:bin_lo).')));
        expCounts(2:nbins-1) = n * (lambda.^(bin_lo+1:bin_hi-1).' * exp(-lambda) ./ factorial((bin_lo+1:bin_hi-1).'));
        expCounts(nbins) = n - sum(expCounts(1:nbins-1));
        if (sum(expCounts < 5) > 0)
%            error('To few counts in a bin.');
        end
        
        % Matlab reference
        if (VECTORIZATION_ON == 1)
            x = poissrnd(lambda,n,1); 
        else
            x = zeros(n,1); for i = 1:n; x(i) = poissrnd(lambda); end;
        end
        
        bins = (bin_lo:bin_hi).';
        obsCounts = zeros(nbins,1); obsCounts(1) = sum(x<=bins(1)); for i = 2:nbins-1; obsCounts(i) = sum(x == bins(i)); end; obsCounts(end) = sum(x>=bins(end));
        if (sum(obsCounts < 5) > 0)
%            error('To few counts in a bin.');
        end
        
        V_ref(k) = sum((obsCounts-expCounts).^2./expCounts);
        
        % Our implementation
        if (VECTORIZATION_ON == 1)
            x = poissrnd_own_parallel(lambda,n,1);
        else
            x = zeros(n,1); for i = 1:n; x(i) = poissrnd_own(lambda); end;
        end
        
        bins = (bin_lo:bin_hi).';
        obsCounts = zeros(nbins,1); obsCounts(1) = sum(x<=bins(1)); for i = 2:nbins-1; obsCounts(i) = sum(x == bins(i)); end; obsCounts(end) = sum(x>=bins(end));
        if (sum(obsCounts < 5) > 0)
%            error('To few counts in a bin.');
        end
        
        V_own(k) = sum((obsCounts-expCounts).^2./expCounts);
    end
    
    if (lambda_idx == 1); 
        figure; 
    end;
    subplot(2,4,lambda_idx); 
    h = cdfplot(V_ref);
    set(h,'Color','b');
    hold on; 
    h = cdfplot(V_own);
    set(h,'Color','k');
    hold on; 
    t = unique(sort([V_ref;V_own]));
    plot(t,chi2cdf(t,dof),'r'); 
    ks_ref = kstest(V_ref,[t,chi2cdf(t,dof)],.05,'unequal');
    ks_own = kstest(V_own,[t,chi2cdf(t,dof)],.05,'unequal');
    title(sprintf('Empirical CDF (\\lambda=%d,\\nu=%d)',lambda,dof));
    legend(sprintf('Matlab (KS=%d)',ks_ref),sprintf('Own (KS=%d)',ks_own),'Theory','Location','SouthEast');
%    [H,P,KSSTAT,CV] = kstest(V_our,[t,chi2cdf(t,dof)],.05,'unequal');
end
