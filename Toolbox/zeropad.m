function Bpad = zeropad(B)
% zeropad with zeros to avoid wrap-around effects in decovolution

[M, N] = size(B);
Bpad = zeros(M+N);
Bpad(int64(M/2)+1:int64(M/2 + M),int64(M/2)+1:int64(M/2 + M)) = B;

end
