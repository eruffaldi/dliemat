function y = se3_fromQuatT(q,t)
q = q/norm(q);
R = quat2mat(q);
det(R)
y =eye(4);
y(1:3,1:3) = R;
y(1:3,4) = t;
end
function R = quat2mat(q)
    sqw = q(1)*q(1);
    sqx = q(2)*q(2);
    sqy = q(3)*q(3);
    sqz = q(4)*q(4);
    
    invs = 1 / (sqx + sqy + sqz + sqw);
    
    m00 = ( sqx - sqy - sqz + sqw)*invs;
    m11 = (-sqx + sqy - sqz + sqw)*invs ;
    m22 = (-sqx - sqy + sqz + sqw)*invs ;
    
    tmp1 = q(2)*q(3);
    tmp2 = q(4)*q(1);
    m10 = 2 * (tmp1 + tmp2)*invs ;
    m01 = 2 * (tmp1 - tmp2)*invs ;
    
    tmp1 = q(2)*q(4);
    tmp2 = q(3)*q(1);
    m20 = 2 * (tmp1 - tmp2)*invs ;
    m02 = 2 * (tmp1 + tmp2)*invs ;
    tmp1 = q(3)*q(4);
    tmp2 = q(2)*q(1);
    m21 = 2 * (tmp1 + tmp2)*invs ;
    m12 = 2 * (tmp1 - tmp2)*invs ;      
    R = [m00 m01 m02;m10 m11 m12;m20 m21 m22];

end