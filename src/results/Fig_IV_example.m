function fig2 = Fig_IV_example(fig2_dir, iX, useFillCases, noPlotCases, pdir, xyLims, LegsPos)
FS_Tick   = 8;
FS_Label  = 9;
Fs_Legend = 7;

% loading data
fig2_file = fullfile(fig2_dir,'processed_data.mat');
load(fig2_file);

% plot figure
[fig2,ax2,axLeg2] = make_fig_m0(m0, iX, useFillCases, noPlotCases,... % plots median & 25th/75th percentile over all seeds
                     LegCols=[3,3],CrameriColors='roma',... -> use subdir w/ only 1 seed per setting for example
                     LegIVnumOnly=true,FigPos=[0.05 0.05 252 300],FigUnits='points',...[-2400 -200 1000 900]
                     FS_Label=FS_Label,FS_Legend=Fs_Legend,LineWidth=1,...
                     Leg_BoxWidth=30); % 15, 12
for k=1:2
    ax2(k).FontSize = FS_Tick;
    ax2(k).YLabel.FontSize = FS_Label;
end
ax2(k).XLabel.FontSize = FS_Label;
fig2.OuterPosition(3:4) = [252 252/fig2.OuterPosition(3)*fig2.OuterPosition(4)];

% getting data from directory
dirs = dir(fig2_dir);
names = {dirs.name};
pattern = sprintf('^0*%d_',iX);
matches = ~cellfun(@isempty, regexp(names, pattern, 'once')); %match pattern
firstMatchIdx = find(matches, 1, 'first'); % Find the first match
pX_dir = fullfile(fig2_dir,names{firstMatchIdx});
cd(pX_dir);
seed_files = dir('seed_*.mat');
cd(pdir);
load(fullfile(pX_dir,seed_files(1).name),'opts','y0','e0','u0'); % makes sense, also only one seed per Re in subdir1

% overlaying e_k & y_k on plots
plot(ax2(1),e0(opts.p+1:end),'r-','DisplayName','$e_k$','LineWidth',1);
plot(ax2(1),y0(opts.p+1:end),'-', 'DisplayName','$y_k$','LineWidth',1,'Color',0.5*ones(1,3));
ax2(1).YLim = xyLims{2};
legax2_1= legend(ax2(1),'Interpreter','latex','FontSize',Fs_Legend,...
                'IconColumnWidth',15,...
                'Position',LegsPos{1});

% overlaying u_k on plots
plot(ax2(2),u0(opts.p+1:end),'-','DisplayName','$u_k \Rightarrow U_{\mathrm{f}}$','LineWidth',1,'Color',0.5*ones(1,3));
ax2(2).YLim = xyLims{3};
legax2_2 = legend(ax2(2),'Interpreter','latex','FontSize',Fs_Legend,...
                'IconColumnWidth',15,...
                'Position',LegsPos{3});

% shifting p steps forwards -> 'future' IO data
offset = opts.p; % shift right by p timesteps
offsetXAxis(ax2(1),offset)
offsetXAxis(ax2(2),offset)

xlim(ax2(1),xyLims{1});

axLeg2(1).Position(1:2) = LegsPos{2};
axLeg2(2).Position(1:2) = LegsPos{4};

end

%% Helper functions

function offsetXAxis(ax, x_offset)
%OFFSETXAXIS Shifts all X data (lines + fills) on an axis by a given offset
%
%   offsetXAxis(gca, 10)   % shifts everything right by +10 units

    % --- Shift Line objects ---
    hLines = findall(ax, 'Type', 'Line');
    for h = hLines'
        h.XData = h.XData + x_offset;
    end

    % --- Shift Patch / Fill objects ---
    hPatches = findall(ax, 'Type', 'Patch');
    for h = hPatches'
        % Case 1: has XData (fill created by 'fill' or 'area')
        if isprop(h, 'XData') && ~isempty(h.XData)
            h.XData = h.XData + x_offset;

        % Case 2: generic patch with Vertices field
        elseif isprop(h, 'Vertices') && ~isempty(h.Vertices)
            verts = h.Vertices;
            verts(:,1) = verts(:,1) + x_offset; % shift x-column
            h.Vertices = verts;
        end
    end

    % Optional: update axis limits
    ax.XLim = ax.XLim + x_offset;
end