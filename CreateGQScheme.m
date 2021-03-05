function [ gq ] = CreateGQScheme(N)
%CreateGQScheme Creates GQ Scheme of order N
%
%   CREATEGQSCHEME creates an array of the Gauss Quadrature scheme of order
%   N for up to fifth order polynomials at N = 3. GQ.npts is the number of
%   Gauss points, GQ.gsw is the weighting of each point, and GQ.xipts is
%   the value of each point.

% Creates and initialises a data structure
gq.npts = N;
if (N > 0) && (N <= 5)
    %order of quadrature scheme i.e. %number of Gauss points
    gq.gsw = zeros(N,1); %array of Gauss weights
    gq.xipts = zeros(N,1); %array of Gauss points
    switch N
        case 1 % first order polynomials
            gq.gsw(1) = 2;
            gq.xipts(1) = 0;
        case 2 % third order polynomials
            gq.gsw(:,1) = 1;
            gq.xipts(1,1) = -sqrt(1/3);
            gq.xipts(2,1) = sqrt(1/3);
        case 3 % fifth order polynomials
            gq.gsw(1,1) = 5/9;
            gq.gsw(2,1) = 8/9;
            gq.gsw(3,1) = 5/9;
            gq.xipts(1,1)=-sqrt(3/5);
            gq.xipts(2,1)=0;
            gq.xipts(3,1)=sqrt(3/5);
        case 4 % seventh order polynomials
            gq.gsw(1,1) = (18+sqrt(30))/36;
            gq.gsw(2,1) = (18+sqrt(30))/36;
            gq.gsw(3,1) = (18-sqrt(30))/36;
            gq.gsw(4,1) = (18-sqrt(30))/36;
            gq.xipts(1,1)=sqrt(3/7-(2/7)*sqrt(6/5));
            gq.xipts(2,1)=-sqrt(3/7-(2/7)*sqrt(6/5));
            gq.xipts(3,1)=sqrt(3/7+(2/7)*sqrt(5/6));
            gq.xipts(4,1)=-sqrt(3/7+(2/7)*sqrt(5/6));
        case 5 % ninth order polynomials
            gq.gsw(1,1) = 128/225;
            gq.gsw(2:3,1) = (322+13*sqrt(70))/900;
            gq.gsw(4:5,1) = (322-13*sqrt(70))/900;
            gq.xipts(1,1)=0;
            gq.xipts(2,1)=(1/3)*sqrt(5-(2)*sqrt(10/7));
            gq.xipts(3,1)=-(1/3)*sqrt(5-(2)*sqrt(10/7));
            gq.xipts(4,1)=-(1/3)*sqrt(5+(2)*sqrt(10/7));
            gq.xipts(5,1)=(1/3)*sqrt(5+(2)*sqrt(10/7));
    end
    
else % higher order polynomials than ninth
    fprintf('Invalid number of Gauss points specified');
end
end