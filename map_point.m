function [ u_,v_,z_ ] = map_point(u, v, h)
    z_ = h(3,1) * u + h(3,2) * v + h(3,3) ;
    u_ = (h(1,1) * u + h(1,2) * v + h(1,3)) ./ z_ ;
    v_ = (h(2,1) * u + h(2,2) * v + h(2,3)) ./ z_ ;
end

