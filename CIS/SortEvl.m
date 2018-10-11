% 
% [evl_r, evl_l] = SortEvl(evl, Nsub, NExtra)
%
% Input:
%   evl(1:Nsub)     -- desired eigenvalues
%   evl(Nsub+1:end) -- extra eigenvalues
%   NExtra          -- number of extra eigenvalues to return (MP: do not need?)
%
% Output:
%   evl_rs  -- evl(1:Nsub) in descending order by real part
%   evl_ls  -- the first Ncompl remaining values (by descending real part)

function [evl_rs, evl_ls] = SortEvl(evl, Nsub)  %MP

evl_r     = evl(1:Nsub);
[tlr, is] = sort(-real(evl_r));
evl_rs    = evl_r(is);

evl_l     = evl(Nsub+1:end);
[tls, is] = sort(-real(evl_l));
evl_ls    = evl_l(is);

