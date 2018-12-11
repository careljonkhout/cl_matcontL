function [q0, q1, p1, p0] = getBTvectors(CISdata, varargin)

global cds

if nargin > 1
    q0 = varargin{1};
    q1 = varargin{2};
    p1 = varargin{3};
    p0 = varargin{4};
else
    
    NSub   = CISdata.NSub;
    Asub   = CISdata.T(1:NSub,1:NSub);
    A0     = CISdata.A;
    Q0     = CISdata.Q(:,1:NSub);
    
    if isequal(cds.curve, @bogdanovtakensL)
        borders = cds.borders; % DV: assures smooth a2 and b2 along BT curve
    else
        [V, D] = eig(Asub);
        [~, ind] = min(abs(diag(D)));
        borders.v = V(:, ind);
        
        [W, D] = eig(Asub');
        [~, ind] = min(abs(diag(D)));
        borders.w = W(:, ind);
    end

    Bord    = [Asub borders.w;borders.v' 0];
    bunit   = [zeros(NSub,1); 1];
    vext    = Bord\bunit;
    v1      = vext(1:end-1);
    q0      = Q0 * v1;
    
    b2 = vext; b2(end) = 0;
    vext2    = Bord\b2;
    v2       = vext2(1:end-1);
    q1       = Q0 * v2;
    
    Bord2    = [A0' Q0*borders.v; borders.w'* Q0' 0];
    bunit2   = [zeros(length(A0),1); 1];
    Wext     = Bord2\bunit2;
    p1       = Wext(1:end-1);
    
    b22 = Wext; b22(end) = 0;
    Wext2    = Bord2\b22;
    p0       = Wext2(1:end-1);
end

% normalize so that <p0,q0>=<p1,q1>=1, <p0,q1>=<p1,q0>=0 and <q0,q0>=1, <q0,q1>=0
mu = sqrt(abs(q0'*q0));
q0 = (1/mu)*q0;
q1 = (1/mu)*q1;
q1 = q1 - (q0'*q1)*q0;
nu = q0'*p0;
p1 = (1/nu)*p1;
p0 = p0 - (p0'*q1)*p1;
p0 = (1/nu)*p0;