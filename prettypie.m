function varargout = prettypie(varargin)

%well, ok, what are the minimum input arguments we need? 
%There are multiple ones. Let's just start with the cell array of arrays,
%since that's closest to what we have
%

% X = avg_P_O3;
% inargs = 
[ax,inargs] = axescheck(varargin);
X = inargs{1}{1};
ax = newplot(ax); 
%%
%First case: X is a cell array of arrays
% X = avg_P_O3;
pie_array = nan(1,100); %Pre-allocate for speed
pie_categories = nan(1,100);
i = 1;
for ind = 1:numel(X)
    curr_array = X{ind}; %Curr array should be a vector, will need to test that
    n = numel(curr_array);
    %Will also need to add some math here for cutoffs
    pie_array(i:i+n-1) = sort(curr_array(:));
    pie_categories(i:i+n-1) = ind;
    i = i + n;
end
pie_array = pie_array(1:i-1);
pie_categories = pie_categories(1:i-1);


%%
%First, make the pie:
h  = pie(pie_array./sum(pie_array)); hold on;

%Ok, now we need to make it look good. 
%We do that by stepping through each pie slice and modifying it. 
starting_cat = 1;
starting_sInd = 1;
% curr_order = all_orders{starting_cat};
all_maps = {'tempo','matter','turbid','speed','amp','-gray'};
curr_map_name = all_maps{starting_cat};
cmap = cmocean(curr_map_name);
l_colors = [];

% Step through once to set the labels
label_cutoff = 1;
for i = 2:2:numel(h)

    if pie_array(i/2)./sum(pie_array) < label_cutoff
            set(h(i),'String','');
        else
            ind = curr_order(starting_sInd);
            set(h(i),'String',components{starting_cat}{ind});


            set(h(i),'fontsize',18);
            pos = get(h(i),'Position');
            angle = atand(pos(2)./pos(1));

            pos;
            r = 1.05;
            new_pos = [abs(cosd(angle)).*pos(1)./abs(pos(1)),abs(sind(angle)).*pos(2)./abs(pos(2)),0].*r;
            set(h(i),'Position',new_pos);


            set(h(i),'Rotation',angle);
            if pos(1) < 0
                set(h(i),'HorizontalAlignment','Right');
            else
                set(h(i),'HorizontalAlignment','Left');
            end
            if i == 6
    %             a(0);
            end


    end 
end

%Step through again to set the colors
%Some math to get the colors right:
for i = 1:max(pie_categories)
    n(i) = sum(pie_categories == i);
end

for i = 1:numel(pie_categories)
    curr_category = pie_categories(i);
    number_within_category = n(curr_category);
    position_within_category = i - sum(n(1:curr_category-1)); %A little weird, but sum(n(1:0)) = 0, so it works
    if position_within_category == 1
        curr_map_name = all_maps{curr_category};
        cmap = cmocean(curr_map_name);
        
        vertices = h(2*i-1).Vertices;
        plot([vertices(end-1,1),vertices(end,1)],[vertices(end-1,2),vertices(end,2)],'-w','linewidth',3);
    end
    
    set(h(2*i-1),'FaceColor',map_colors(cmap,[0,number_within_category+1],position_within_category));
    set(h(2*i-1),'EdgeColor',map_colors(cmap,[0,number_within_category+1],position_within_category));
end
figure(h(1).Parent.Parent);


%%
a = 18;