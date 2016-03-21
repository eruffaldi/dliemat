% SE3 parametric form Gaussian
%
% TODO: add derivative with gaussian
% TODO: add integral with gaussian
%
% Emanuele Ruffaldi 2015
classdef se3d
    properties
        x %compacted as for the matse3 se3d_*
    end
    
    methods
        function this = se3d(mu,Sigma)
            if isa(mu,'se3d') & nargin == 1
                this.x = mu.x;
            elseif length(mu) == (36+16)
                this.x = mu(:)'; % flattened
            else
                assert(ndims(mu) == 2 & diall(size(mu) == 4),'wrong mean size expected 4x4');
                if nargin == 1
                    Sigma = zeros(6); % exact
                else
                    assert(ndims(mu) == 2 & all(size(Sigma) == 6),'wrong covariance size expected 6x6');
                end
                this.x = se3d_set(mu,Sigma);
            end
        end
        
        function n = length(this)
            n = 6;
        end
        
        function m = mean(this)
            [m,S] = se3d_get(this.x);
        end
        
        function this = interp(this,other,alpha)
            this.x = se3d_interp(this.x,other.x,alpha);
        end

        % variance
        function S = var(this)
            [m,S] = se3d_get(this.x);
        end

        % std
        function C = std(this)
            [m,S] = se3d_get(this.x);
            C = chol(S);
        end
        
        % inversion
        function this = inv(this)
            this.x = se3d_inv(this.x);
        end
    
        % logarithm of the mean
        % TODO: implement a class for the uncertain velocity: se3dvel
        function r = log(this)
            [m,S] = se3d_get(this.x);
            r = se3_log(m);
        end

        % minus means distance
        function d = minus(this,other)
            d = se3_dist(se3d_mean(this.x),se3d_mean(other.x));
        end

        % -x = inverse
        function this = uminus(this,other)
            this.x = se3d_inv(this.x);
        end

        % se3ggauss fuse se3ggauss = se3ggauss
        function this = fuse(this,other)
            this.x = s3d_fuse(this.x,other.x);
        end

        %samples the distribution n times => returns n copies of the group
        function r = sample(this,n)
            r = se3d_sample(this.x,n);
        end
        
        % se3d * se3d => se3d 
        % se3d * groupvar(4x4) => se3d
        % se3d * flattengauss  => se3d
        % groupvar(4x4) * se3d => se3d TODO
        function this = mtimes(this,other)
            if isa(this,'se3d')
                if isa(other,class(this))
                    this.x = se3d_mul(this.x,other.x);    
                else
                    if all(size(other) == [4,4])
                        this.x = se3d_mul(this.x,se3d_set(other,zeros(6)));
                    elseif length(other) == (16+36)
                        this.x = se3d_mul(this.x,other);
                    end
                end
            else
                error('other * se3d not implemented')
            end
        end

    end
    
end

