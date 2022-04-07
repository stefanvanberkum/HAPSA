# HAPSA
Health-adjusted assortment planning and shelf-space allocation (HAPSA).

## Set-up instructions.
- Set-up CPLEX for Java:
  - Make sure you have a recent version of Java JDK installed on your device (we used Java 15.0.2), if not download from: https://www.oracle.com/nl/java/technologies/javase-downloads.html (most recent version) or https://www.oracle.com/java/technologies/javase/jdk15-archive-downloads.html (Java 15.0.2).
  - You can open and edit the code in any editor, but we recommend the Eclipse Java IDE (we used 2021-03): https://www.eclipse.org/downloads/.
  - Download and install IBM ILOG CPLEX Optimization Studio: https://www.ibm.com/academic/technology/data-science -> Software.
  - Download the software from this repository and open in Eclipse: File -> Open Projects from File System...
  - Set-up CPLEX for the Eclipse Java IDE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-eclipse-java-api.

- Tune CPLEX parameters (optional):
  - When running the algorithm, it is recommended to use tuned CPLEX parameters: either tuned to your specific problem or the ones pre-tuned by us.
  - For tuning: run ParameterTuner.java for your specific problem (an example of how we tuned our data is provided in the main).
  - Please note that tuning might take a very long time depending on your problem and problem size. If so, consider using the pre-tuned parameters or a time limit for each tuning run.

## How to use?
- Run Main.java for your specific problem (the way we ran it is provided in the main).
  - Any store can be simulated using the StoreSimulator.java class.
  - The possible objective types are: Assortment Planning and Shelf-space Allocation (APSA), APSA with an availability penalty (AVA), Health-adjusted APSA (HAPSA), APSA with a Healthy-Left, Unhealthy Right approach (HLUR), and APSA with a visibility penalty (VIS).
- For each model, a summary is written to a text file and the decision variables s_kj, x_ij, and y_kj are written to a CSV file.
- Results can be plotted using plotResults.R (we ran it in RStudio).
- Shelf colormaps can be obtained by running plotShelf.m (we ran it in MATLAB R2020a).

## How to cite?
BibTeX:

@thesis{VanBerkum2021, 
author = {{Van Berkum}, Stefan}, 
school = {Erasmus University Rotterdam}, 
title = {{Factoring health into the equation: Promoting healthy purchasing decisions through shelf-space optimization}}, 
type = {bachelor's thesis}, 
url = {http://hdl.handle.net/2105/59932}, 
year = {2021}}

APA:

Van Berkum, S. (2021). _Factoring health into the equation: Promoting healthy purchasing decisions through shelf-space optimization._ (bachelor's thesis). Erasmus University Rotterdam. Retrieved from http://hdl.handle.net/2105/59932
