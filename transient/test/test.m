function test()
handle1 = @diff;
handle2 = @add;
save handle handle1 handle2
cell = {handle1 handle2};
feval(cell(1),4,3)
end
function [sum] = add(a,b)
sum = a + b;
sum
end

function [dif] = diff(a,b)
dif = a - b;
dif
end

