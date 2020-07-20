%%% This function generates random variables from a triangular distribution.
%%% Input: a = min, b = mode, c = max, n = number of RVs
%%% Output: X - column random vector
function X = TriRnd(a,b,c,n)
U = rand(n,1);

% use Inverse Transform Method
% F^-1(u) = a + sqrt((b-a)*(c-a))*sqrt(u),   if u in [0, (b-a)/(c-a)]
%         = c - sqrt((c-b)*(c-a))*sqrt(1-u), otherwise
k1 = sqrt((b-a)*(c-a));
k2 = sqrt((c-b)*(c-a));
part1 = a + k1 * sqrt(U);
part2 = c - k2 * sqrt(1-U);

I = U <= (b-a)/(c-a);
X = part1 .* I + part2 .* (~I);

% % test example
% a = 0; b = 3; c = 5; n = 100000;
% X = TriRnd(a,b,c,n);
% [mean(X) (a+b+c)/3]  % compare with the true mean
% [var(X) (a^2 + b^2 + c^2 - a*b - a*c - b*c)/18] % compare with the true variance