%% 标准前向列推法
function [P]=legendre_col_zh(degree,theta)
%theta=90-lambda;
p=zeros(degree+1);
p(1,1)=1;
p(2,1)=sqrt(3)*cosd(theta);
p(2,2)=sqrt(3)*sind(theta);
[r,c]=size(p);
u=sind(theta);
t=cosd(theta);

%n=m
for i=3:r
    n=i-1;
    p(i,i)=u*sqrt((2*n+1)/(2*n))*p(i-1,i-1);
    %当m=n-1时，bnm=0
    p(i,i-1)=sqrt(2*n+1)*t*p(i-1,i-1);
end

%n>m
for i=3:r
    n=i-1;
    for j=1:i-2
        m=j-1;
        anm=sqrt(((2*n-1)*(2*n+1))/((n-m)*(n+m)));
        bnm=sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((n-m)*(n+m)*(2*n-3)));
        p(i,j)=anm*t*p(i-1,j)-bnm*p(i-2,j);
    end
end

calall=(1+degree)^2/2;
P=zeros(calall,3);
d=1;
for i=1:r
    for j=1:i
        P(d,1)=i-1;
        P(d,2)=j-1;
        P(d,3)=p(i,j);
        d=d+1;
    end
end
end

