# categorypie
CATEGORYPIE Displays a pie chart organized by category

CATEGORYPIE(X,G) draws a pie chart of the data in the vectory X grouped by
the categorical array G

CATEGORYPIE(C), where C is a cell array of arrays, draws a pie chart of all
the data in C grouped by each array in C [e.g. C = {[1,2,3], [3,3,2],
[1,4,1.5]} ]

CATEGORYPIE(S), where S is a structure, draws a pie chart of all the data
in S grouped by the data in S [e.g. S = struct('first_group,[1,2,3],
'second_group',[3,3,2], 'third_group', [1,4,1.5]) ] 

CATEGORYPIE(AX, ...) produces a pie chart in axes with handle AX.

CATEGORYPIE(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
parameter name/value pairs.

Optional input variables (passed as parameter/value pairs): [default]
  plotcutoff:         Scalar. Minimum fraction of the total for a slice to be 
                      included in the pie chart. Any entries less than 
                      this fraction are removed before plotting.  [0]

  labelcutoff:        Scalar. Minimum fraction of the total for a slice to be 
                      labeled. Has no effect if slices are not labeled.
                      [0.01]

  categorylabels:     Cell Array. Labels for each category of pie slices. Required if
                      labelmode is 'category', unless the inputs were given
                      as as structure. Must be the same size as the number 
                      of categories. []

  slicelabels:        Cell Array of Cells or Cell Array or Structure of cells; 
                      Labels for each individual pie slice. Required if
                      labelmode is 'slice'. The form and size of
                      slicelabels must match that of the input data.
                      []. 

  sorted:             Boolean. If true, the order of slices in each
                      category will be plotted from smallest to largest.
                      [true]

  rotatelabels:       Boolean. If true, labels will be rotated to be
                      aligned  with the slice they are labeling. Useful
                      to labeling small slices [true]. 

  labelfontsize:      Scalar. Font size to plot the labels. [Matches the
                      font size of the current axes]

  labelmode:          String. Specifies how the pie should be labeled.
                      Must be: 'none', 'auto','category','slice',
                      'percentage'. ['auto']
                          'none': No labels
                          'percentage': Percentage of total
                          'category': One label for each category of slices
                          'slice': One label for each slice
                          'auto': Chooses based on contexts. Picks the first of 
                                 'slice' -> 'category' ->  'percentage' that is
                                  possible.

  usecmocean:         Boolean. If true, will use cmocean colormaps
                      (Unless overridden by 'colormaps'), if false will
                      use MATLAB built in colormaps. 

  colormaps:          Cell Array. Lists the colormaps to color each 
                      category. Each entry in the cell array can be
                      a valid entry to colormap.m, or a valid entry to
                      cmocean.m. Named entries to colormap.m can be
                      prepended with a dash to reverse them, e.g.
                      '-parula'. [{'tempo', 'matter', 'turbid', 'speed',
                       'amp','-gray', '-ice', '-pink', 'algae'} if cmocean is available,
                      otherwise {'-summer','-copper','-bone','-pink','-gray','-hot'}]

  colormap#:          Any valid entry to colormap.m or cmocean.m. '#' can
                      be any number, and the entry given here will
                      replace that entry in colormaps. 




Output variables:
  h:                  Vector of text and patch handles, in alternating
                      order. Equivalent to the output from pie.m

  category_patches:   Vector of patches, one to represent each category.
                      Convenient for representing each category in a
                      legend. 

  all_maps:           Cell array of the colormaps (each [nx3]) used in
                      to color each category. 
  
  edges:              Vector of patch handles to the white dividing lines
                      between categories. 



Example Usages:

Cell-Based input:
  categorypie({[1,2,3],[3,3,2],[1,4,1.5]},'labelcutoff',0,'labelmode','auto','categorylabels',{'Ordered','Random','Decimal'});

Array Based Input:
  categorypie([1,2,3,3,3,2,1,4,1.5],categorical([1,1,1,2,2,2,3,3,3]),
      'slicelabels',{'A1','A2','A3','B1','B2','B3','C1','C2','C3'},'sorted',false,'usecmocean',false)
  
Structure Based Input:
  categorypie(struct('Ordered',[1,2,3],'Random',[3,3,2],'Decimal',[1,4,1.5]))
