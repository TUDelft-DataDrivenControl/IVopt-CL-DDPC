function fig = Fig_Lf_estimates(Cases,mLf,iX,opts,notScaledto)
FS_Tick  = 8;
FS_Label = 9;

nCases = numel(Cases);
iXstr = sprintf('iX%d',iX);

Lf_act = mLf.actLf.(iXstr);
maxabsval = 0;
for k = 1:nCases
    switch Cases{k}
        case 'actLf'
            Lf_plt = Lf_act;
        case notScaledto
            continue; % do not scale clim to account for these IVs
        otherwise
            Lf_plt = mLf.(Cases{k}).(iXstr).mean - Lf_act;
    end
    maxabsval = max(max(abs(Lf_plt),[],'all'),maxabsval);
end

fig = figure();
axs = gobjects(nCases,1);
fig.OuterPosition(4) = ceil(108.8333*nCases);

% Adjust margins for subplots to allow space for separator bars
space4line   = 0.04;
space_axes   = 0.02;
space_bottom = 0.03;
subplot_height = (1-space4line-space_bottom-(nCases-2)*space_axes)/nCases;

% calculate y positions for subplots
yPositions = calc_yPositions(nCases,subplot_height,space4line,space_axes,space_bottom);

for k = 1:nCases
    switch Cases{k}
        case 'actLf'
            Lf_plt = Lf_act;
            yLabel = 'actual';
        otherwise
            Lf_plt = mLf.(Cases{k}).(iXstr).mean - Lf_act;
            if strcmp(Cases{k},'TrPred')
                yLabel = 'TP';
            elseif strcmp(Cases{k},'CLSPC')
                yLabel = 'CL-SPC';
            else
                yLabel = ['IV',Cases{k}(3:end)];
            end            
    end

    % create subplot in determined position
    axs(k) = subplot('Position',[0.1, yPositions(k), 0.85, subplot_height]);

    imagesc(Lf_plt);
    clim([-maxabsval maxabsval]);
    cmap = crameri('vik','pivot',0);
    colormap(cmap);
    xline(opts.p*opts.nu+0.5); xline(opts.p*(opts.nu+opts.ny)+0.5)
    ylabel(yLabel,'Interpreter','latex','FontSize',FS_Label);
    axs(k).FontSize = FS_Tick;
    if k~=nCases
        axs(k).XTickLabel = {};
    end
end

% Colorbar
fig.Units = 'points';
fig.OuterPosition(3:4) = [252 252/fig.OuterPosition(3)*fig.OuterPosition(4)];
fig.Units = 'normalized';
drawnow;

% align bottom axes
yPositions = align_bottom(axs,yPositions);

% update axes position based on bottom
TopPos    = sum(axs(1).OuterPosition([2, 4]));
BottomPos = axs(end).OuterPosition(2);
yStrechFac = 1/(TopPos - BottomPos);
subplot_height = subplot_height * yStrechFac;
space4line     = space4line     * yStrechFac;
space_axes     = space_axes     * yStrechFac;
space_bottom   = space_bottom   * yStrechFac;
yPositions = calc_yPositions(nCases,subplot_height,space4line,space_axes,space_bottom);
for k = 1:nCases
    axs(k).Position(2) = yPositions(k);
end
drawnow;

% align bottom axes
yPositions = align_bottom(axs,yPositions);

% calculate position of new colorbar
xCB_spacing = 0.02; % spacing between sublots and colorbar
xCB = sum(axs(end).Position([1,3])) + xCB_spacing;
yCB = axs(end).Position(2);
wCB = 1 - xCB;
hCB = sum(axs(1).Position([2,4]))-axs(end).Position(2);

% Add colorbar manually to the right of the last subplot
cb = colorbar(axs(end),'Location','manual','FontSize',FS_Tick,'Position',[xCB, yCB, wCB, hCB]);
cb.Location = 'manual';
cb.Position = [xCB, yCB, wCB, hCB];%axCB.Position;
drawnow;

% room left over on the right
pos = getCBpos(cb);
RoomRight = 1-sum(pos([1,3]));
for ii = 1:10
    if RoomRight ~= 0
        RightTotalRoom = 1 - axs(1).Position(1);
        RoomNeeded = RightTotalRoom - RoomRight;
        xScaleFac = RightTotalRoom/RoomNeeded;
        for k = 1:nCases
            axs(k).Position(3) = axs(k).Position(3)*xScaleFac;
        end
        xCB_spacing = xCB_spacing*xScaleFac;
        xCB = sum(axs(end).Position([1,3])) + xCB_spacing;
        wCB = wCB*xScaleFac;
        cb.Position = [xCB, yCB, wCB, hCB];
    end
    drawnow;
    pos = getCBpos(cb);
    RoomRight = 1-sum(pos([1,3]));
    if abs(RoomRight) <= 1e-4
        break;
    end
end

end

%% Local functions
function cbOuterPos = getCBpos(cb)

    % Create temporary axes to measure colorbar position
    originalUnits = cb.Units;
    measAxes = axes('Position', cb.Position, 'YAxisLocation', 'right', ...
            'YLim', cb.Limits, 'Units', cb.Units, 'YTick', cb.Ticks,...
            'FontWeight', cb.FontWeight, 'Visible', 'on', ...
            'FontSize', cb.FontSize, 'XTick', [], 'FontName', cb.FontName, ...
            'TickDir',cb.TickDirection, 'YTickLabels', cb.TickLabels);

    if ~isempty(cb.Label.String)
        ylabel(measAxes, cb.Label.String, 'FontSize', cb.Label.FontSize, ...
        'FontWeight', cb.Label.FontWeight)
    end

    % Adjust for (estimated) additional space between colorbar and tick labels
    measAxes.Units = 'pixels';
    measAxes.Position(3) = measAxes.Position(3) + 4; % add 4 pixels (adjust as needed)
    measAxes.Units = originalUnits;
    
    cbOuterPos = get(measAxes, 'OuterPosition');
    delete(measAxes);
end

% calculate y positions for subplots
function yPositions = calc_yPositions(nCases,subplot_height,space4line,space_axes,space_bottom)
    yPositions = nan(nCases,1);
    bottom = 1;
    for k2 = 1:nCases
        if k2 == 1
            spacer = space4line;
        elseif k2 == nCases
            spacer = space_bottom;
        else
            spacer = space_axes;
        end
        yPositions(k2) =  bottom - subplot_height;
        bottom = yPositions(k2) - spacer;
    end
end

% make bottom axes fit
function yPositions = align_bottom(axs,yPositions)
    nCases = numel(axs);
    BottomPos = axs(end).OuterPosition(2);
    axs(end).OuterPosition(2) = 0;
    for k = 1:nCases-1
        axs(k).Position(2) = axs(k).Position(2) - BottomPos;
    end
    drawnow;
    yPositions = yPositions - BottomPos*ones(size(yPositions));
end