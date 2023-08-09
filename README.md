# RShiny App for Bioinformatic Analysis 
This application is a user-friendly RShiny application designed to simplify bioinformatic analyses of mRNA sequencing data. Tailored for researchers, biologists, and bioinformaticians, this Bioinformatics Analysis webapp aims to streamline the process of understanding genetic information and extracting meaningful insights.

## Key Features 
1. **Effortless Data Upload and Visualisation**: Simply upload the files from your computer and visualise the data in a searchable table in  the application.
2. **Sample Information Exploration**: Load and examine sample information matrix. View distribution of all continous variables in the sample information matrix as density plots.
3. **Counts Matrix Exploration**: Input normalized counts matrix by any method as a CSV and be able to choose different gene filtering thresholds and assess their effects using diagnostic plots of the counts matrix. Available analyses include:
   - Diagnostic scatter plots of median count vs variance and median count vs number of zeros.
   - Clustered heatmap of counts remaining after filtering.
   - Scatter plot of principal component analysis projections where user can visualize the top *N* principal components in a beeswarm plot.
4. **Explore Differential Expression**: Simply upload differential expression dataset to visualise it as a volcano plot with customizable X axis and Y axis. View the filtered DEGs in a sortable table.
5. **Visualization of Individual Gene Expression(s)**: Selects counts from an arbitrary gene and visualize them broken out by a desired sample information variable as a bar plot, boxplot, violin plot, or beeswarm plot.

## Future Directions
1. Develop a Test Suite to test all back-end functions.
2. Allow user to upload raw data and perform normalization and differential expression analysis in the back-end for a more user friendly experience. 
