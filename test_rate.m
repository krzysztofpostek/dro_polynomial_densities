for deg = 1:4  
    for r = 1:20
        %[bound] = Lasserre_heuristic_v2(r,1);
        [bound] = Lasserre_bound_box(r,deg)
        value(deg,r) = (bound +1)*(r^2);
    end

end


for deg = 1:4
    plot(1:20,value(deg,:),'-+')
    hold on
end    