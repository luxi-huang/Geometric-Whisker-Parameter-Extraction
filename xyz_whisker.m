function [x,y,z] = xyz_whisker(S,A,N)


s = linspace(0,S,N);
x = zeros(1,length(s));
y = zeros(1,length(s));
z = zeros(1,length(s));
for j = 1:N
    x(j) = x_from_a_and_s(A,s(j));
    z(j) = -A*x(j)^2;
end