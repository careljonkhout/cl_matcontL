function [vv,dd] = cdf2rdf(v,d)
i = find(imag(diag(d))');
index = i(1:2:length(i));
if isempty(index)
   vv=v; dd=d;
else   
   if (max(index)==size(d,1)) | any(conj(d(index,index))~=d(index+1,index+1))
      error(message('MATLAB:cdf2rdf:invalidDiagonal'));
   end
   j = sqrt(-1);
   t = eye(length(d));
   twobytwo = [1 1;j -j];
   for i=index
      t(i:i+1,i:i+1) = twobytwo;
   end 
   vv=v/t; dd=t*d/t;
end
