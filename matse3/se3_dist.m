function d = se3_dist(g1,g2)

d = se3_log(se3_mul(se3_inv(g1),g2));