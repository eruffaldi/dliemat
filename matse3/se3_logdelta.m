function omega = se3_logdelta(g1,g2)

omega = se3_log(se3_mul(g1,se3_inv(g2)));