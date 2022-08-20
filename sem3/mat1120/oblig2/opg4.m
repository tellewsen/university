function [ C ] = opg4( a,b )
n = length(a);
if length(a) == length(b)
    for i=1:n
        for j=1:n
            A(i,j) = a(i)^(j-1);
        end
    end
else
    disp('vektor a og b må ha lik lengde')
    return
end
if det(A) == 0
    disp('A ble ikke inverterbar')
    return
else
    C = inv(A)*b;
end
end
%Kjøreeksempel:
%opg4([1;1],[1;1])
%A ble ikke inverterbar
%opg4([1;1;1],[1;1])
%vektor a og b må ha lik lengde
%opg4([1;1],[1;1;1])
%vektor a og b må ha lik lengde
%opg4([1;5],[3;7])
%
%ans =
%
%     2
%     1