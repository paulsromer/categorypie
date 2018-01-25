function varargout = categorypie(varargin)
%CATEGORYPIE Displays a pie chart organized by category
%
% CATEGORYPIE(X,G) draws a pie chart of the data in the vectory X grouped by
% the categorical array G
%
% CATEGORYPIE(C), where C is a cell array of arrays, draws a pie chart of all
% the data in C grouped by each array in C [e.g. C = {[1,2,3], [3,3,2],
% [1,4,1.5]} ]
%
% CATEGORYPIE(S), where S is a structure, draws a pie chart of all the data
% in S grouped by the data in S [e.g. S = struct('first_group,[1,2,3],
% 'second_group',[3,3,2], 'third_group', [1,4,1.5]) ] 
%
% CATEGORYPIE(AX, ...) produces a pie chart in axes with handle AX.
%
% CATEGORYPIE(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
% parameter name/value pairs.

% Optional input variables (passed as parameter/value pairs): [default]
%   plotcutoff:         Scalar. Minimum fraction of the total for a slice to be 
%                       included in the pie chart. Any entries less than 
%                       this fraction are removed before plotting.  [0]
% 
%   labelcutoff:        Scalar. Minimum fraction of the total for a slice to be 
%                       labeled. Has no effect if slices are not labeled.
%                       [0.01]
% 
%   categorylabels:     Cell Array. Labels for each category of pie slices. Required if
%                       labelmode is 'category', unless the inputs were given
%                       as as structure. Must be the same size as the number 
%                       of categories. []
%
%   slicelabels:        Cell Array of Cells or Cell Array or Structure of cells; 
%                       Labels for each individual pie slice. Required if
%                       labelmode is 'slice'. The form and size of
%                       slicelabels must match that of the input data.
%                       []. 
%
%   sorted:             Boolean. If true, the order of slices in each
%                       category will be plotted from smallest to largest.
%                       [true]
%
%   rotatelabels:       Boolean. If true, labels will be rotated to be
%                       aligned  with the slice they are labeling. Useful
%                       to labeling small slices [true]. 
%
%   labelfontsize:      Scalar. Font size to plot the labels. [Matches the
%                       font size of the current axes]
%
%   labelmode:          String. Specifies how the pie should be labeled.
%                       Must be: 'none', 'auto','category','slice',
%                       'percentage'. ['auto']
%                           'none': No labels
%                           'percentage': Percentage of total
%                           'category': One label for each category of slices
%                           'slice': One label for each slice
%                           'auto': Chooses based on contexts. Picks the first of 
%                                  'slice' -> 'category' ->  'percentage' that is
%                                   possible.
%
%   usecmocean:         Boolean. If true, will use cmocean colormaps
%                       (Unless overridden by 'colormaps'), if false will
%                       use MATLAB built in colormaps. cmocean colormaps
%                       are available at https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps
%
%   colormaps:          Cell Array. Lists the colormaps to color each 
%                       category. Each entry in the cell array can be
%                       a valid entry to colormap.m, or a valid entry to
%                       cmocean.m. Named entries to colormap.m can be
%                       prepended with a dash to reverse them, e.g.
%                       '-parula'. [{'tempo', 'matter', 'turbid', 'speed',
%                        'amp','-gray', '-ice', '-pink', 'algae'} if cmocean is available,
%                       otherwise {'-summer','-copper','-bone','-pink','-gray','-hot'}]
%
%   colormap#:          Any valid entry to colormap.m or cmocean.m. '#' can
%                       be any number, and the entry given here will
%                       replace that entry in colormaps. 
%
%
%
%
% Output variables:
%   h:                  Vector of text and patch handles, in alternating
%                       order. Equivalent to the output from pie.m
%
%   category_patches:   Vector of patches, one to represent each category.
%                       Convenient for representing each category in a
%                       legend. 
%
%   all_maps:           Cell array of the colormaps (each [nx3]) used in
%                       to color each category. 
%   
%   edges:              Vector of patch handles to the white dividing lines
%                       between categories. 
%
%
%
% Example Usages:
%
% Cell-Based input:
%   categorypie({[1,2,3],[3,3,2],[1,4,1.5]},'labelcutoff',0,'labelmode','auto','categorylabels',{'Ordered','Random','Decimal'});
%
% Array Based Input:
%   categorypie([1,2,3,3,3,2,1,4,1.5],categorical([1,1,1,2,2,2,3,3,3]),
%       'slicelabels',{'A1','A2','A3','B1','B2','B3','C1','C2','C3'},'sorted',false,'usecmocean',false)
%   
% Structure Based Input:
%   categorypie(struct('Ordered',[1,2,3],'Random',[3,3,2],'Decimal',[1,4,1.5]))


%%
hg2flag = ~verLessThan('matlab', '8.4.0');
r2016aflag = ~verLessThan('matlab', '9.0.0');
r2013bflag = ~verLessThan('matlab', '8.2.0');
cmoceanflag = exist('cmocean.m','file') == 2;

%-------------------
% Parse Inputs
%-------------------

%Pop out the axes, if any. 
[ax,inargs] = axescheck(varargin);
input_data = inargs{1}{1};
ax = newplot(ax); 

%Then figure out which set of inputs we have here
pv = inargs{1};
if iscell(pv{1}) %First option: Cell array of strings
    cellfun(@(x) validateattributes(x,{'numeric'},{'vector'},1),pv{1});
    nReqArgs = 1;
    input_type = 'cell';
    
elseif isnumeric(pv{1}) %Second option: Numeric Array
    validateattributes(pv{1},{'numeric'},{'nonnegative','real','finite'})
    %If the first input is numeric, then the second input must be
    %categorical
    if ~iscategorical(pv{2})
        error('Numeric Array input to categorypie must be accompanied by a categorical array')
    else
        validateattributes(pv{2},{'categorical'},{'size',size(pv{1})},'categorypie','G'); %size of the categorical array must equal size of the input_data
    end    
    nReqArgs = 2;
    input_dictionary = pv{2};
    input_type = 'numeric';
elseif isstruct(pv{1}) %Third option: Structure organized by category
    validateattributes(pv{1},{'struct'},{'scalar'})
    structfun(@(x) validateattributes(x,{'numeric'},{'vector'}),pv{1}); 
    nReqArgs = 1;
    input_type = 'struct_category';
    
else
    error('Input type not recognized. categorypie.m accepts as input a cell array, a structure, or a numerical array')
end


optArgs = pv(nReqArgs+1:end);

%Now parse the additional inputs. 
if r2013bflag
    addParamMethod = 'addParameter';
else
    addParamMethod = 'addParamValue';
end

p = inputParser;
p.(addParamMethod)('categorylabels',    {},     @(x) validateattributes(x, {'cell'},{}));
p.(addParamMethod)('slicelabels',       {},     @(x) validateattributes(x, {'cell'},{}));
p.(addParamMethod)('labelcutoff',       0.01,   @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative'}));
p.(addParamMethod)('plotcutoff',        0,      @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative'}));
p.(addParamMethod)('sorted',            true,   @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));
p.(addParamMethod)('labelfontsize',     ax.FontSize,    @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive','finite','real'}));
p.(addParamMethod)('labelmode',         'auto', @(x) validateattributes(x, {'string','char'},{}));
p.(addParamMethod)('colormaps',         {},     @(x) validateattributes(x, {'cell'},{}))
p.(addParamMethod)('usecmocean',        [],  @(x) validateattributes(x, {'logical','numeric'},{'scalar'}));
p.(addParamMethod)('rotatelabels',      true,   @(x) validateattributes(x,{'logical','numeric'},{'scalar'}));

p.KeepUnmatched = true;
p.PartialMatching = true;
p.parse(optArgs{:});
Opt = p.Results;
Ex = p.Unmatched;

%Extra logic to validate and parse 'labelmode', 'usecmocean', and
%'colormaps'
validatestring(Opt.labelmode, {'none', 'auto','category','slice','percentage'}, 'categorypie', 'labelmode');

if isempty(Opt.usecmocean) %Don't override use-specified input
    Opt.usecmocean = cmoceanflag;
end
if Opt.usecmocean && ~cmoceanflag
    error('usecmocean flag specified but cmocean.m not found');
end
if isempty(Opt.colormaps) %Don't override user-specified input
    if Opt.usecmocean
        Opt.colormaps = {'tempo','matter','turbid','speed','amp','-gray','-ice','-pink','algae'};
    else
        Opt.colormaps = {'-summer','-copper','-bone','-pink','-gray','-hot'};
    end
end
   


%%
%--------------------------------------------------
% Translate the data into a form pie.m can read
%--------------------------------------------------

pie_array = nan(1,100); 
pie_categories = nan(1,100);
pie_slicelabels = cell(1,100);
use_slicelabels = ~isempty(Opt.slicelabels);


i = 1;
switch input_type 
    case 'cell'
        if use_slicelabels %Check that form of the labels matches the form of the data
            validateattributes(Opt.slicelabels,{'cell'},{'numel',numel(input_data),'vector'})
        end
        for ind = 1:numel(input_data) %input data is a cell array of arrays
            curr_array = input_data{ind};
            n = numel(curr_array);
            [pie_array, I] = append_array(pie_array, curr_array, i, n, Opt);
            pie_categories(i:i+n-1) = ind;
            if use_slicelabels
                curr_labels = Opt.slicelabels{ind};
                assert(numel(curr_labels) == numel(curr_array),sprintf('Number of labels does not match number of entries in category %i',ind))                    
                pie_slicelabels(i:i+n-1) = curr_labels(I);
            end
            i = i + n;
        end
        numer_categories = numel(input_data);
        
        
    case 'numeric'  
        pie_categorylabels = input_dictionary;
        fn = categories(pie_categorylabels);
        if use_slicelabels
            validateattributes(Opt.slicelabels,{'cell'},{'size',size(input_dictionary)});
        end
        for ind = 1:numel(fn)
            m = pie_categorylabels == fn{ind};
            curr_array = input_data(m);
            n = numel(curr_array);
            [pie_array, I] = append_array(pie_array, curr_array, i, n, Opt);
            pie_categories(i:i+n-1) = ind;
            if use_slicelabels
             curr_labels = Opt.slicelabels(m);
             assert(numel(curr_labels) == numel(curr_array),sprintf('Number of labels does not match number of entries in category %i',ind))                    
             pie_slicelabels(i:i+n-1) = curr_labels(I);
            end
            i = i + n;
        end
        numer_categories = numel(fn);
        
    case 'struct_category'
        pie_categorylabels = fieldnames(input_data);
        if use_slicelabels
            validateattributes(Opt.slicelabels,{'struct'},{'scalar'});
            names_match = all(cellfun(@(x1) strcmp(x1,x2),fieldnames(Opt.slicelabels),fieldnames(pie_categorylabels)));
            assert(names_match,'Input Strucuture and slicelabels must have the same fieldnames');
        end
        for ind = 1:numel(pie_categorylabels)
            curr_fieldname = pie_categorylabels{ind};
            curr_array = input_data.(curr_fieldname);
            n = numel(curr_array);
            [pie_array, I] = append_array(pie_array, curr_array, i, n, Opt);
            pie_categories(i:i+n-1) = ind;
            if use_slicelabels
                curr_labels = Opt.slicelabels.(curr_fieldname);
                pie_slicelabels(i:i+n-1) = curr_labels(I);
            end
            i = i+n;
        end
        numer_categories = numel(pie_categorylabels); 
end

m = pie_array./nansum(pie_array) > Opt.plotcutoff; %Apply the plot cutoff math
pie_array = pie_array(m);
pie_categories = pie_categories(m);
if use_slicelabels
    pie_slicelabels = pie_slicelabels(m);
end
numer_entries = numel(pie_array);


%Validate categorylabels 
if ~isempty(Opt.categorylabels)
    validateattributes(Opt.categorylabels,{'cell'},{'numel',numer_categories});
else %If categorylabels is empty, we default to using the categories supplied by the categorical array or the structure
    switch input_type 
        case 'numeric'
            Opt.categorylabels = categories(input_dictionary);
        case 'struct_category'
            Opt.categorylabels = fieldnames(input_data);
    end
end
%After this point, how the inputs are received should not matter

%%
%--------------------------------------------------------------
% Build the pie chart and create the full list of colormaps
%--------------------------------------------------------------
h  = pie(pie_array./sum(pie_array)); hold on;


%% Decide on the label mode and step through once to set the labels
if strcmp(Opt.labelmode,'auto')
    if ~isempty(Opt.slicelabels)
        Opt.labelmode = 'slice';
    else
        if ~isempty(Opt.categorylabels)
            Opt.labelmode = 'category';
        else
            Opt.laleblmode = 'percentage';
        end
    end
end

for i = 2:2:numel(h)
    set(h(i),'fontsize',Opt.labelfontsize);
    switch Opt.labelmode
        case 'none'
            set(h(i),'String','');
        case 'percentage'
            if pie_array(i/2)./sum(pie_array) < Opt.labelcutoff
                set(h(i),'String','');
            end
        case 'category'
            if isempty(Opt.categorylabels)
                error('categorypie.m called with labelmode option ''category'' but no category labels given');
            end
            set(h(i),'String',''); %There's more to this section later on. 
    
            
        case 'slice'
            if isempty(Opt.slicelabels)
                error('categorypie.m called with labelmode option ''slice'' but no slice labels given');
            end
            if pie_array(i/2)./sum(pie_array) < Opt.labelcutoff
                set(h(i),'String','');
            else
                [slice] = adjust_label(h(i), pie_slicelabels{i/2},Opt.rotatelabels)
%                 set(h(i),'String',pie_slicelabels{i/2});
%                 
%                 pos = get(h(i),'Position');
%                 angle = atand(pos(2)./pos(1));
% 
%                 pos;
%                 r = 1.05;
%                 new_pos = [abs(cosd(angle)).*pos(1)./abs(pos(1)),abs(sind(angle)).*pos(2)./abs(pos(2)),0].*r;
%                 set(h(i),'Position',new_pos);
% 
% 
%                 set(h(i),'Rotation',angle);
%                 if pos(1) < 0
%                     set(h(i),'HorizontalAlignment','Right');
%                 else
%                     set(h(i),'HorizontalAlignment','Left');
%                 end
            end
    end
end
%%
%---------------------------------
% Assign colors to each slice
%---------------------------------

all_maps = Opt.colormaps;
if numer_categories > numel(all_maps)
    all_maps = repmat(all_maps,1,ceil(numer_categories./numel(all_maps)));
    all_maps = all_maps(1:numer_categories);
end

fn = fieldnames(Ex); %Now try to find any individually specified colormaps
for ind = 1:numel(fn)
    curr_fn = fn{ind};
    if strncmp('colormap',curr_fn,8)
        number = str2double(curr_fn(9:end));
        name = Ex.(curr_fn);
        all_maps{number} = name;
    else
        warning(['Optional parameter ', curr_fn, ' not recognized.']);
    end
end

%Some math to get the colors right:
for i = 1:max(pie_categories)
    n(i) = sum(pie_categories == i);
end

edges = [];
for i = 1:numel(pie_categories)
    curr_category = pie_categories(i);
    number_in_category = n(curr_category);
    position_within_category = i - sum(n(1:curr_category-1));
    if position_within_category == 1
        curr_map_name = all_maps{curr_category};
        cmap = assign_cmap(curr_map_name);
        
        %Create a white boundary between categories
        %Make the boundaries a patch to make them adjust as the figure size
        %changes
        vertices = h(2*i-1).Vertices;
        theta1 = atan2(vertices(2,2),vertices(2,1));
        theta2 = theta1 + .015;
        theta = [theta1:1e-3:theta2, theta2];
        X = [0, cos(theta), 0]; Y = [0, sin(theta), 0];
        edges(i) = patch(X,Y,'w','facecolor','w','edgecolor','w','linewidth',2);  
    end
    set(h(2*i-1),'FaceColor',map_colors(cmap,[0,number_in_category+1],position_within_category));
    set(h(2*i-1),'EdgeColor',map_colors(cmap,[0,number_in_category+1],position_within_category));
end


%% 
%-----------------------------------------------------------------------------
% Pick a slice to act as the representative for each category and label if
% necessary
%-----------------------------------------------------------------------------

%Select slices that are late in the category because they have darker
%colors
for i = 1:max(pie_categories)
    curr_category = i;
    number_in_category = n(curr_category);
    shift = sum(n(1:curr_category-1));
    if number_in_category <= 3 
        %For categories with few slices, pick the last slice
        category_representative_ind(i) = number_in_category + shift;
    else
        %Selects the largest slice that's the late-middle of the
        %cateogry
        ind = floor(0.6*number_in_category):floor(0.9*number_in_category);
        [c, indind] = max(pie_array(ind+shift));
        category_representative_ind(i) = ind(indind)+shift;
    end
end

category_patches = 2*category_representative_ind-1; 
%Use those labels if the labelmode is 'category'
if strcmp(Opt.labelmode,'category')
    for ind = 1:numel(category_representative_ind)
        m = pie_categories == ind;
        total_of_that_category = sum(pie_array(m))./sum(pie_array);
        curr_ind = 2*category_representative_ind(ind);
        if total_of_that_category > Opt.labelcutoff
            [slice] = adjust_label(h(curr_ind), Opt.categorylabels{ind},Opt.rotatelabels)
        end
    end
end

%%
%----------------------------
% Assign output variables
%----------------------------
figure(h(1).Parent.Parent);

out = {h, category_patches, all_maps, edges};
vargout = out(1:nargout);

%%
%-------------------------
% Subfunctions
%-------------------------
function [new_array, I] = append_array(old_array,curr_array, i, n, Opt)
    if Opt.sorted
        [old_array(i:i+n-1), I] = sort(curr_array(:));
    else
        old_array(i:i+n-1) = curr_array(:);
        I = 1:n;
    end
    new_array = old_array;
    
function [cmap] = assign_cmap(cmap_name)
    try
        cmap = cmocean(cmap_name);
    catch E
        if strcmp(E.message,'Unrecognized colormap name.') || strcmp(E.identifier,'MATLAB:UndefinedFunction')
            if strncmp(cmap_name,'-',1)
                cmap_name = cmap_name(2:end);
                cmap = colormap(cmap_name);
                cmap = flipud(cmap);
            else
                cmap = colormap(cmap_name);
            end
        else
            rethrow(E)
        end
    end
    
    
function  [slice] = adjust_label(slice, new_text, rotate_labels)
    set(slice,'String',new_text);

    pos = get(slice,'Position');
    angle = atand(pos(2)./pos(1));

    pos;
    r = 1.05;
    new_pos = [abs(cosd(angle)).*pos(1)./abs(pos(1)),abs(sind(angle)).*pos(2)./abs(pos(2)),0].*r;
    set(slice,'Position',new_pos);

    if rotate_labels
        set(slice,'Rotation',angle);
    end
    if pos(1) < 0
        set(slice,'HorizontalAlignment','Right');
    else
        set(slice,'HorizontalAlignment','Left');
    end
    
function colors = map_colors(cmap,c_ax,input_vals)
    n = size(cmap,1);
    c_range = linspace(c_ax(1),c_ax(2),n);

    trimmed_vals = input_vals;
    trimmed_vals(trimmed_vals < c_ax(1)) = c_ax(1);
    trimmed_vals(trimmed_vals > c_ax(2)) = c_ax(2);
    colors = nan(numel(trimmed_vals),3);
    for ind = 1:numel(trimmed_vals)
        curr_tv = trimmed_vals(ind);
        ci = find(curr_tv >= c_range,1,'last');
        colors(ind,:) = cmap(ci,:);
    end
    
       
    
    
        
