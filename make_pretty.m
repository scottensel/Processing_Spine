function make_pretty(f)
% make_pretty does basic figure clean-up
%
% It can be called directly after plotting
% e.g. 
%     >>  figure; hold on;
%     >>  plot(randn(10,2),randn(10,2),'.'); 
%     >>  make_pretty; 
%
% INPUTS
%     f (optional) : figure handle
%
% make_pretty.m adjusts a few properties
%     1. Makes figure larger and sets background color to 'w'
%     2. Sets axes tick directions to 'out', makes them larger, and gray, and
%        increases font size
%     3. Sets marker size ('.' to 15, 'o' and others to 5); 
%     4. Sets line width to 1
%     5. Sets EdgeColor to 'none' (for bar or patch) and FaceAlpha to 0.5
%        (for e.g. patch)
%     6. Adjusts legend if it exists
%%
if nargin==0 || isempty(f)
    f = gcf; % get figure handle
end

% 1: Resize 
if all(f.Position(3:4)==[560 420])
    f.Position(1:2) = f.Position(1:2).*.8; 
    f.Position(3:4) = [560 420]*1.25; 
end
set(f,'Color','w'); 

subplot_num = length(f.Children); % Number of axes/subplots
for q = 1:subplot_num % loop through axes
    if contains(char(class(f.Children(q))),'Axes')
        h = f.Children(q); % axis handle

        % 2: Fix Ticks and Axes
        set(h,'TickDir','out'); 
        box(h,'off'); 
        set(h,'TickLength',[.02 .02]); 
        set(h,'FontSize',14); 
        set(h,'LabelFontSizeMultiplier',1.25); 
        set(h,'XColor',[0.25 0.25 0.25]); 
        set(h,'YColor',[0.25 0.25 0.25]); 

        % 3: Set MarkerSize
%         for i = 1:length(h.Children)
%             if isprop(h.Children(i),'MarkerSize')
%                 if strcmp(h.Children(i).Marker,'.')
%                     h.Children(i).MarkerSize = 15; 
%                 else
%                     h.Children(i).MarkerSize = 5; 
%                 end
%             end
%         end
        % 4: Set LineWidth
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'LineWidth') && h.Children(i).LineWidth==0.5
                h.Children(i).LineWidth = 1; 
            end
        end

        % 5: Set EdgeColor and FaceAlpha
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'EdgeColor') 
                has_markerface = isprop(h.Children(i),'MarkerFaceColor') && ~strcmp(h.Children(i).MarkerFaceColor,'none'); 
                has_face = isprop(h.Children(i),'FaceColor') && ~(strcmp(h.Children(i).FaceColor,'none') || strcmp(char(h.Children(i).FaceColor),char([0 0 0]))); 
                if has_markerface||has_face
                    h.Children(i).EdgeColor = 'none'; 
                end
            end
        end
        for i = 1:length(h.Children)
            if isprop(h.Children(i),'FaceAlpha') 
                if strcmp(get(h.Children(i),'Type'),'bar')
                    h.Children(i).FaceAlpha = 0.9;
                else
                    h.Children(i).FaceAlpha = 0.5;
                end
            end
        end

        % 6: Check if legend exists and clean it up
        if ~isempty(h.Legend)
            leg = h.Legend; 
            set(leg,'Location','best'); 
            set(leg,'Box','off'); 
            set(leg,'FontSize',14); 
        end
    end
end