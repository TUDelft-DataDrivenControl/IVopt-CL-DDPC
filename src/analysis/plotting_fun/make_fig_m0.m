function make_fig_m0(m0, pX, opts)
%MAKE_FIG_M0  Prepare and plot IV trajectories with bounds for a given index.
%
%   make_fig_m0(m0, pX, opts)
%
%   Inputs:
%     m0   Struct with fields Uf and Yf, each containing IV data.
%     pX   Index for the iX string (integer).
%     opts Options struct passed to plot_IV_trajectories.

iX_str = sprintf('iX%d',pX);
% select the 25th and 75th percentile -> 6 & 16th of 21
% see size(m0.Uf.iv1.iX1.pctiles,3);

Uf_ivs = fieldnames(m0.Uf);
Uf_lb_fs = replace(Uf_ivs,'iv','l');
Uf_ub_fs = replace(Uf_ivs,'iv','u');
Uf_av_fs = replace(Uf_ivs,'iv','m');
num_Uf_ivs = numel(Uf_ivs);

u_iv = struct;
for kIV = 1:num_Uf_ivs
    lb_name = Uf_lb_fs{kIV};
    ub_name = Uf_ub_fs{kIV};
    av_name = Uf_av_fs{kIV};
    iv_name = Uf_ivs{kIV};
    
    u_iv.(av_name) = m0.Uf.(iv_name).(iX_str).median;
    u_iv.(lb_name) = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,6));
    u_iv.(ub_name) = squeeze(m0.Uf.(iv_name).(iX_str).pctiles(:,:,16));
end

Yf_ivs = fieldnames(m0.Yf);
Yf_lb_fs = replace(Yf_ivs,'iv','l');
Yf_ub_fs = replace(Yf_ivs,'iv','u');
Yf_av_fs = replace(Yf_ivs,'iv','m');
num_Yf_ivs = numel(Yf_ivs);

y_iv = struct;
for kIV = 1:num_Yf_ivs
    lb_name = Yf_lb_fs{kIV};
    ub_name = Yf_ub_fs{kIV};
    av_name = Yf_av_fs{kIV};
    iv_name = Yf_ivs{kIV};
    
    y_iv.(av_name) = m0.Yf.(iv_name).(iX_str).median;
    y_iv.(lb_name) = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,6));
    y_iv.(ub_name) = squeeze(m0.Yf.(iv_name).(iX_str).pctiles(:,:,16));
end
plot_IV_trajectories(u_iv,y_iv,opts,"useBounds",true);
end
