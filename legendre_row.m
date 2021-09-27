%% 标准前向行推法
function [P]=legendre_row_zh(degree,theta)
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
end

%n>m
for i=3:r
    n=i-1;
    for k=i-1:-1:1
        m=k-1;
        if m==0
            j=2;
        else
            j=1;
        end
        gnm=(2*(m+1))/sqrt((n+m+1)*(n-m));
        hnm=sqrt(((n+m+2)*(n-m-1))/((n-m)*(n+m+1)));
        % m=n-1,hnm=0
        if hnm==0
            p(i,k)=(1/sqrt(j))*(gnm*t/u*p(i,k+1));
        else
            p(i,k)=(1/sqrt(j))*(gnm*t/u*p(i,k+1)-hnm*p(i,k+2));
        end
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