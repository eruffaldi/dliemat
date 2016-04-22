%
% Sets the Group in SE3 from position and quaternion xyzw
%
% Emanuele Ruffaldi SSSA 2015-2016
function y = se3_setpq(p,q_xyzw)

y = se3_set(quat_to_rot(q_xyzw),p);

function R = quat_to_rot( Qrotation )
% qGetR: get a 3x3 rotation matrix
% R = qGetR( Qrotation )
% IN: 
%     Qrotation - quaternion describing rotation
% 
% OUT:
%     R - rotation matrix 
%     
% VERSION: 03.03.2012


x = Qrotation( 1 );
y = Qrotation( 2 );
z = Qrotation( 3 );
w = Qrotation( 4 );

Rxx = 1 - 2*(y^2 + z^2);
Rxy = 2*(x*y - z*w);
Rxz = 2*(x*z + y*w);

Ryx = 2*(x*y + z*w);
Ryy = 1 - 2*(x^2 + z^2);
Ryz = 2*(y*z - x*w );

Rzx = 2*(x*z - y*w );
Rzy = 2*(y*z + x*w );
Rzz = 1 - 2 *(x^2 + y^2);

R = [ 
    Rxx,    Rxy,    Rxz;
    Ryx,    Ryy,    Ryz;
    Rzx,    Rzy,    Rzz];