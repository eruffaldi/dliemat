function y = se3_fromQuatT(q,t)
q = q/norm(q);
R = quat2mat(q);

y =eye(4);
y(1:3,1:3) = R;
y(1:3,4) = t;
end

function [a,theta] = quat2axis(q)

a = q(2:4);
theta = 2*atan2(norm(a),q(1));
a = a / sin(theta); % == sqrt(1-w*w)
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

% 
% %
% % Sets the Group in SE3 from position and quaternion xyzw
% %
% % Emanuele Ruffaldi SSSA 2015-2016
% function y = se3_setpq(p,q_xyzw)
% 
% y = se3_set(quat_to_rot(q_xyzw),p);
% 
% function R = quat_to_rot( Qrotation )
% % qGetR: get a 3x3 rotation matrix
% % R = qGetR( Qrotation )
% % IN: 
% %     Qrotation - quaternion describing rotation
% % 
% % OUT:
% %     R - rotation matrix 
% %     
% % VERSION: 03.03.2012
% 
% 
% x = Qrotation( 1 );
% y = Qrotation( 2 );
% z = Qrotation( 3 );
% w = Qrotation( 4 );
% 
% Rxx = 1 - 2*(y^2 + z^2);
% Rxy = 2*(x*y - z*w);
% Rxz = 2*(x*z + y*w);
% 
% Ryx = 2*(x*y + z*w);
% Ryy = 1 - 2*(x^2 + z^2);
% Ryz = 2*(y*z - x*w );
% 
% Rzx = 2*(x*z - y*w );
% Rzy = 2*(y*z + x*w );
% Rzz = 1 - 2 *(x^2 + y^2);
% 
% R = [ 
%     Rxx,    Rxy,    Rxz;
%     Ryx,    Ryy,    Ryz;
%     Rzx,    Rzy,    Rzz];