function trun_Angles = truncate_Angle(Angles, a, b, type_left, type_right)
% =========================================================================
% -- Function to truncate angle values in a specific range (a, b) depending
%    on the "open" or "close" sides of the interval;
%    this function supports input 'Angles' as multi-dimensional arrays 
% =========================================================================

% -- Function: Angles_trun = truncate_Angle(Angles, a, b, type_left, type_right)

% -- Input Arguments:
%       Angles: Input angles, a 3D-Array of size(p.nClusters, p.nRays, num_freqs_per_subcarr)
%       a: Minimum "left" side of the interval
%       b: Maximum "right" side of the interval
%       type_left: Options are: 'open', 'close'
%       type_right: Options are: 'open', 'close'

% -- Output Arguments:
%       trun_Angles: Truncated angles based on the input angles and the interval limits,
%        a 3D-Array of size(p.nClusters, p.nRays, num_freqs_per_subcarr).

%=================================================

% -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, 
%             Mohamed Habib Loukil, Hakim Jemaa, 
%             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

% -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
%            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
%            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa

% =========================================================================

% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO:  A  channel  simulator for  wideband  ultra-massive  MIMO  terahertz  communications," 
%  arXivpreprint arXiv:2104.11054, 2021.

% =========================================================================

trun_Angles = Angles;
EPSILON = 1e-5;

% left side
if strcmp(type_left,'open')
    for indx_clusters = 1: size(Angles,1)       % p.nClusters
        for indx_rays = 1: size(Angles,2)       % p.nRays
            for indx_freq = 1: size(Angles,3)   % p.nFreq(1) = num_freqs_per_subcarr
                
                if(Angles(indx_clusters,indx_rays,indx_freq) <= a)
                    trun_Angles(indx_clusters,indx_rays,indx_freq) = a+EPSILON;
                end
                
            end
        end
    end
else %close
    for indx_clusters = 1: size(Angles,1)
        for indx_rays = 1: size(Angles,2)
            for indx_freq = 1: size(Angles,3)
                
                if(Angles(indx_clusters,indx_rays,indx_freq) < a)
                    trun_Angles(indx_clusters,indx_rays,indx_freq) = a;
                end
                
            end
        end
    end
end

% right side
if strcmp(type_right,'open')
    for indx_clusters = 1: size(Angles,1)
        for indx_rays = 1: size(Angles,2)
            for indx_freq = 1: size(Angles,3)
                
                if(Angles(indx_clusters,indx_rays,indx_freq) >= b)
                    trun_Angles(indx_clusters,indx_rays,indx_freq) = b-EPSILON;
                end
                
            end
        end
    end
else %close
    for indx_clusters = 1: size(Angles,1)
        for indx_rays = 1: size(Angles,2)
            for indx_freq = 1: size(Angles,3)
                
                if(Angles(indx_clusters,indx_rays,indx_freq) > b)
                    trun_Angles(indx_clusters,indx_rays,indx_freq) = b;
                end
                
            end
        end
    end
end

end