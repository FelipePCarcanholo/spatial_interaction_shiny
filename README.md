# spatial_interaction_shiny
shiny app to visualize spatial data with scores interacted

## Using my data to test shiny
To use the app with the data I provided, download this 2 dataframes in this drive:
https://drive.google.com/drive/folders/1fPeKKG-RUh3WDDyim14D83Nt7QzvQ6Na?usp=sharing

Confirm that is installed all the packages needed in the begining of .R file
Finaly run all the code in .R file

## Using your own data in the app
You will use only metadata, not need to use seurat object

Make sure that the spatial location of x and y are in the metadata with the columns written as "x" and "y"

#### You must:
1. Put metadata on the argument "metadata"
2. Put the name of the column which has the slice/patient information in the argument "patient_column", using "" in the name
3. Put a vector of names of the columns that contain the scores that you will analyse in the argument "list_of_scores_columns"
4. Put the name of the column with the cell type information (or cluster if you want) in the argument "cell_type_column", using "" in the name
5. Put the colors for each cell type in object "color" in the begining of code

#### You can:
1. Put the name of the column that separate the type of tumor between patients (to separate the plot in the final part that uses all patients to make the median score of interaction), its in the argument "proximity_split", using "" in the name. If you leave this part with just "" it will plot for all patients in the metadata.
2. Select how many neighbohrs its supposed to look in the median interaction spatial score, defoult its 5.


## App functionality
### 1. Filter data
Fist you will select which slice you want to analyse, can only analyse one per time

### 2. Select the scores you want to plot
The names of columns in "list_of_scores_coluns" will appear as options of check box, you select the ones you want to analyse (can be how many you want)

### 3. Click in Generate Plots & Heatmaps 
This will generate:

1. A Spatial plot colored by the cell type you had specified
This plot its interactive, so you can select an area of interest to plot the heatmap of all scores and the mean spatial interaction score for just the cells you have selected in that area

2. A Spatial plot colored by the intensity of the score (one for each score you have selected in the checked boxes)

3. A Heatmap of the mean score of the cells interacted (one for each score you have selected in the checked boxes), written inside of the squares the pvalue of that score (based on 1k simulations of resampling randomly scores for that interaction), pvalue <0.05 is colored by green and >0.05 in black

4. A Heatmap of all scores for all cells for that slice you selected (separated by cell type)

### 4. You can select a range of rows (scores) from the overall Heatmap 
This will plot a filtred version of the overall Heatmap 
(usefull when you have too many scores to analyse in the overall heatmap)

### 5. You can click on "All Samples Score Interaction"
This will plot the mean score of spatially interacted cells considering all patients, separated by the column you have put on "proximity_split", to obtain the mean score of interaction for that specific tumor you are analysing. Written inside of the squares the pvalue of that score (based on 1k simulations of resampling randomly scores for that interaction), pvalue <0.05 is colored by green and >0.05 in black
