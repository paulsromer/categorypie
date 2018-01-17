function varargout = prettypie(varargin)

%well, ok, what are the minimum input arguments we need? 
%There are multiple ones. Let's just start with the cell array of arrays,
%since that's closest to what we have
%

hg2flag = ~verLessThan('matlab', '8.4.0');
r2016aflag = ~verLessThan('matlab', '9.0.0');
r2013bflag = ~verLessThan('matlab', '8.2.0');

%-------------------
% Parse Inputs
%-------------------

%Pop out the axes, if any. 
[ax,inargs] = axescheck(varargin);
input_data = inargs{1}{1};
ax = newplot(ax); 

%Then figure out which set of inputs we have here
pv = inargs{1};
if iscell(pv{1})
    %First option: We enter a cell array of strings
    %Requirements: All entries in the cell array must be one-dimensional
    %numeric, all real and positive

    %Using validateattributes for this
    cellfun(@(x) validateattributes(x,{'numeric'},{'vector'},1),pv{1}) %,'nonnegative','real','finite','nonnan'
    nReqArgs = 1;
    input_type = 'cell';
    
elseif isnumeric(pv{1})
    %If the first input is numeric, then it must satisfy all the important
    %constraints:
    validateattributes(pv{1},{'numeric'},{'nonnegative','real','finite'})
    
    %If the first input is numeric, then the second input must be
    %categorical
    if ~iscategorical(pv{2})
        error('Numeric Array input to prettypie must be accompanied by a categorical array')
    else
        validateattributes(pv{2},{'categorical'},{'size',size(pv{1})}); %size of the categorical array must equal size of the input_data
    end    
    nReqArgs = 2;
    input_dictionary = pv{2};
    input_type = 'numeric';
elseif isstruct(pv{1}) %Require struct by category
    validateattributes(pv{1},{'struct'},{'scalar'})
    structfun(@(x) validateattributes(x,{'numeric'},{'vector'}),pv{1}); %,'nonnegative','real','finite','nonnan'
    nReqArgs = 1;
    input_type = 'struct_category';
    
        
else
    error('Input type not recognized. prettypie.m accepts as input a cell array, a structure, or a numerical array')
end


optArgs = pv(nReqArgs+1:end);

%Now parse the additional inputs. 
if r2013bflag
    addParamMethod = 'addParameter';
else
    addParamMethod = 'addParamValue';
end

%what if we ditched the labels, and instead had 'categorylabels' and
%'slicelabels'? %And then also had a 'labelmethod' and a 'legendmethod'
%options? That might work better. 

p = inputParser;
% p.addOptional('labels',             {},     @(x) validateattributes(x, {'cell'}, {}));
p.(addParamMethod)('categorylabels',    {}, @(x) validateattributes(x, {'cell'},{}));
p.(addParamMethod)('slicelabels',       {}, @(x) validateattributes(x, {'cell'},{}));
p.(addParamMethod)('labelcutoff',   1,      @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative'}));
p.(addParamMethod)('plotcutoff',    0,      @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative'}));
p.(addParamMethod)('sorted',        true,  @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
p.(addParamMethod)('labelfontsize',  ax.FontSize,    @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive','finite','real'}));
p.(addParamMethod)('colorscheme',  {},      @(x) validateattributes(x, {'cell'}, {}));
p.(addParamMethod)('labelmode',    'auto',  @(x) validateattributes(x, {'string','char'},{}));
p.(addParamMethod)('trimnans',      false,  @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));

p.KeepUnmatched = true;
p.parse(optArgs{:});
Opt = p.Results;
a = 18;

validatestring(Opt.labelmode, {'none', 'auto','category','slice','percentage'}, 'prettypie', 'labelmode');
%%
%Options for labels:
%if input_type == 'cell' -> then slicelabels should also be a cell of cells with the same
%size; categorylabels should be a cell
%if input_type == 'array' -> then slicelabels should be the same size as
%array
%if input_type == 'structcat' -> then slicelabels should be the same size
%as structcat
%So I think I want to remove structIndiv. 
%%
pie_array = nan(1,100); %Pre-allocate for speed
pie_categories = nan(1,100);
pie_slicelabels = cell(1,100);
        
%First case: input_data is a cell array of arrays
switch input_type
    case 'cell'
       
        i = 1;
        for ind = 1:numel(input_data)
            curr_array = input_data{ind};
            n = numel(curr_array);
            if Opt.sorted
                [pie_array(i:i+n-1), I] = sort(curr_array(:));
            else
                pie_array(i:i+n-1) = curr_array(:);
                I = 1:n;
            end
            pie_categories(i:i+n-1) = ind;

            curr_labels = Opt.slicelabels{ind};
            pie_slicelabels(i:i+n-1) = curr_labels(I);
            i = i + n;
        end
        m = pie_array./nansum(pie_array) > Opt.plotcutoff;
        pie_array = pie_array(m);
        pie_categories = pie_categories(m);
        pie_slicelabels = pie_slicelabels(m);
        pie_categorylabels = Opt.categorylabels;
    case 'numeric'
        pie_categorylabels = input_dictionary;
        fn = categories(pie_categorylabels);
        i = 1;
        for ind = 1:numel(fn)
            m = pie_categorylabels == fn{ind};
            curr_array = input_data(m);
            n = numel(curr_array);
            if Opt.sorted
                [pie_array(i:i+n-1), I] = sort(curr_array(:));
            else
                pie_array(i:i+n-1) = curr_array(:);
                I = 1:n;
            end
            pie_categories(i:i+n-1) = ind;
            i = i + n;
%             curr_labels = Opt.slicelabels(m);
%             pie_slicelabels(i:i+n-1) = ind;
        end
        m = pie_array./nansum(pie_array) > Opt.plotcutoff;
        pie_array = pie_array(m);
        pie_categories = pie_categories(m);
%         pie_slicelabels = pie_slicelabels(m);
        pie_categorylabels = Opt.categorylabels;        
    case 'struct_category'
        pie_categorylabels = fieldnames(input_data);
        i = 1;
        for ind = 1:numel(pie_categorylabels)
            curr_fieldname = pie_categorylabels{ind};
            curr_array = input_data.(curr_fieldname)
            n = numel(curr_array);
            if Opt.sorted
                [pie_array(i:i+n-1), I] = sort(curr_array(:));
            else
                pie_array(i:i+n-1) = curr_array(:);
                I = 1:n;
            end
            pie_categories(i:i+n-1) = ind;
            i = i+n;
        end
        m = pie_array./nansum(pie_array) > Opt.plotcutoff;
        pie_array = pie_array(m);
        pie_categories = pie_categories(m);
        pie_categorylabels = pie_categorylabels;
end

%After this point, it shouldn't matter how the inputs were recieved, and we
%can focus on displaying the pie correctly. 
% if trim

%%
%First, make the pie:
h  = pie(pie_array./sum(pie_array)); hold on;

%Ok, now we need to make it look good. 
%We do that by stepping through each pie slice and modifying it. 
starting_cat = 1;
starting_sInd = 1;
all_maps = {'tempo','matter','turbid','speed','amp','-gray'};
curr_map_name = all_maps{starting_cat};
cmap = cmocean(curr_map_name);
l_colors = [];

% Step through once to set the labels
if strcmp(Opt.labelmode,'auto')
    %Some logic here to figure out the label mode
end

for i = 2:2:numel(h)
    if pie_array(i/2)./sum(pie_array) < Opt.labelcutoff
        set(h(i),'String','');
    else
%         ind = curr_order(starting_sInd);
%         set(h(i),'String',components{starting_cat}{ind});
        set(h(i),'String',pie_slicelabels{i/2});

        set(h(i),'fontsize',Opt.labelfontsize);
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

a = 18;
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
%         plot([vertices(1,1),vertices(2,1)],[vertices(1,2),vertices(2,2)],'-w','linewidth',3);
        theta1 = atan2(vertices(2,2),vertices(2,1));
        theta2 = theta1 + .015;
        theta = [theta1:1e-3:theta2, theta2];
        X = [0, cos(theta), 0]; Y = [0, sin(theta), 0];
        patch(X,Y,'w','facecolor','w','edgecolor','w','linewidth',2);
        
    end
    
    set(h(2*i-1),'FaceColor',map_colors(cmap,[0,number_within_category+1],position_within_category));
    set(h(2*i-1),'EdgeColor',map_colors(cmap,[0,number_within_category+1],position_within_category));
end
figure(h(1).Parent.Parent);


%%
a = 18;