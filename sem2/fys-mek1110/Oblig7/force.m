function F = force (dx ,d,k,q)
if dx <d
F = k* abs (dx -d).^q;
else
F = 0.0;
end