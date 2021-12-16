Locust Trajectory Data via Video Tracking
===

This repository contains code used in creating the work:

* J Weinburd, J Landsberg, A Kravtsova, S Lam, T Sharma, SJ Simpson, GA Sword, J Buhl. *Anisotropic interaction and motion state of locusts in a hopper band.*
Preprint on bioRxiv.
<https://doi.org/10.1101/2021.10.29.466390>

The code here collects and analyzes numerical position-time data from raw video footage of locust hoppers in the field. The video was recorded by Drs. Stephen J Simpson, Greg A Sword, and Jerome Buhl and provided as part of an ongoing collaboration. Students who have contributed to this project include: Anna Kravtsova (Eastern Washington U), Shanni Lam (HMC), **Jacob Landsberg** (Haverford), and Tarush Sharma (HMC).

There are three levels of functionality:
1. Download only this repository and immediately run scripts in `\Figures` and `\Tables` using the stored data contained therein.
2. Download the full data set <https://doi.org/10.5061/dryad.n02v6wwzz>, allowing you to recreate the figures and tables directly from the final data set, and to run the script `pipeline_examples.m` to clean and process a few example data sets and compute the tracking accuracy.
3. Install Fiji on your own computer and set up the following (details and links in References at the end). This allows you to run the example pipeline to extract trajectory data from the videos provided with the data downloaded in #2.
    * the TrackMate v6.0.3 plugin
    * the Matlab-ImageJ update site

There are a few files are in the main directory:

* `README.md`
* `pipeline_examples.m`
    * represents a full implementation of the data extraction, cleaning, and processing
    * calls routines in \Functions in sequence, starting from either:
        * videos (requires Fiji with TrackMate v6.0.3), or 
        * existing TrackMate data files
    * cleans and processes the manual tracking data via same routines
    * assesses the accuracy of the automatic tracking
    * saves processed data to data_recording_examples.mat
* `video_preprocess.py`
    * calls `ffmpeg` commands to clip and deinterlace the original video file
    * for exposition only, since the full video files are not publicly available
* `supp_vid2.m`
    * overlays the processed tracking data on top of one of the 10sec example clips and saves a .mpeg video.

The remainder of this README describes the sub-directories and their contents.

\Data
---

Download the project's data from:
<https://doi.org/10.5061/dryad.n02v6wwzz>

Unpack the file into this directory. The directory should now contain Matlab data files `.mat` which include cleaned and processed numerical trajectories and two subdirectories.

1. `data_recording.mat`, `data_recording2.mat`
    * contains the the full data set used in the paper

**\examples** contains all data files for the two example clips analyzed in `pipeline_examples.m`. Includes the starting preprocessed video, intermediate processed video and `.xml` data files, and final data `data_recording_examples.mat`.

**\manual** contains the two ground-truth data sets obtained by manual tracking in TrackMate. **Jacob Landsberg** and **Tarush Sharma** did the bulk of the manual tracking. The original video files are also included so that the `.xml`s can be opened in Fiji>Plugins>Tracking>Load a TrackMate File.


\Figures
---

Each script loads `fig_data.mat` or reassembles the data necessary from the full data set `data_recording.mat` and `data_recording2.mat`.

1. `fig_data.mat` includes data required for all figures.
2. `fig2.m`
    * plots density, polarization, and entropy as in Figure 2
3. `fig3.m`
    * plots histograms of instantaneous speed as appear in Figure 3 and Figure 1 (left).
4. `fig4.m`
    * plots relative neighbor density as appears in Figure 4 and Figure 1 (right).
5. `fig5.m`
    * plots measures of anisotropy as a function of distance from the focal individual as appears in Figure 5.

\Tables
---

Each script loads `table_data.mat` or reassembles the data necessary from the full data set `data_recording.mat` and `data_recording2.mat`.

1. `table_data.mat` includes data required for all tables in the appendix.
2. `app_table1.m`
    * builds a Matlab table with values for density, polarization, and entropy as in Appendix Table 1.
3. ``app_table2.m`
    * builds a Matlab table with percent locusts in each motion state and measures of speed as in Appendix Table 2.
4. ``app_table3.m`
    * builds a Matlab table with measures of anisotropy in the distribution of nearest neighbors as in Appendix Table 3.


\Functions
---

Contains Matlab functions that are called in the process of extracting, cleaning, and processing data. All inputs/outputs in typewriter font are Matlab variables used in `pipeline_examples.m`.

1. data2struct.m
    * **input:** `data_final`, `data_neighbors`
    * **output:** `data_struct`
    * we convert our data to struct arrays before we save it to reduce size
    * Matlab inconveniently forces you to compress .mat files > 2 GB, uncompressing to load a 2.7 GB variable actually takes a while
2. dataClean.m
    * relies on a custom version of a script that ships with TrackMate `trackMateEdges_JWB.m`
    * **input:** `data_rough`
    * optional input to treat data as coming from manual tracking; that is, to include motion state
    * **output:** `data_clean`
    * rewritten from a script written by **Jacob Landsberg** (Haverford College)
3. dataForm.m
    * **input:** `data_clean`,`model` (a trained SVM)
    * **output:** `data_final`
    * modified from a script written by **Jacob Landsberg** (Haverford College)
4. dataNeighbors.m
    * **input:** `data_final`
    * **ouput:** `data_neighbors`
    * computes the relative position of all neighbors for all locusts
    * uses a particle-in-cell method for speed
5. fitSVM.m
    * **input:** none, but does load `data_recording_examples.mat`
    * **output:** a fitcecoc model
    * trains a SVM and on Jacob's manually tracked data
6. initData.m
    * is called by pipeline.m
    * no **input** or **output**
    * loads recording.mat and reads video metadata into it
7. linkAccuracy.m
    * **input1:** `data_clean`, `data_clean_manual`
    * **input2:** `LSC` = Link Similarity Coefficient
8. projTrans.m
    * **input:** coordinates of the corners of the board (in pixels) and the dimensions of the field of view (in pixels)
    * **output:** a 3x3 matrix representing a projective transformation, the average scaling of that transformation on video field of view, and fieldDims = [width, height] of the smallest rectangle that fits the transformed field of view
9. spotAccuracy.m
    * **input1:** `data_clean`, `data_clean_manual`
    * **input2:** `SSC` = Spot Similarity Coefficient
    * modified from a script written by **Jacob Landsberg** (Haverford College)
10. struct2data.m
    * **input:** `data_struct`
    * **output:** `data_final`, `data_neighbors`
11. track.m
    * relies on ImageJ (Fiji), the TrackMate v6.0.3 plugin, and the Matlab-ImageJ update site
    * **input:** video.avi file(s)
    * **output:** a string of the video name, but executes: ImageJ video processing, TrackMate tracking, and saves the result as an .xml file
12. trackmateEdges_JWB.m
    * **input:** data.xml file
    * **output:** a Matlab table of TrackMate edges
    * modified from a script that ships with TrackMate, written by **Jean-Yvez Tinevez**
13. xml2mat.m
    * relies on Matlab scripts that ship with TrackMate
    * **input:** data.xml file
    * optional input to treat data as coming from manual tracking
    * **output:** `data_rough` Matlab variable

**\packages** contains two packages required by functions above. **\TrackMateScripts** includes routines that ship with the TrackMate plugin. **\CircStat** includes a few functions from a package for circular statistics


\Dependencies and References
---

1. ffmpeg
    * FFmpeg Developers. *ffmpeg*. Open Source, Accessed May 2020. <https://ffmpeg.org/ffmpeg.html#Authors>.
2. Fiji (ImageJ2 distribution)
    * Johannes Schindelin, Ignacio Arganda-Carreras, Erwin Frise, Verena Kaynig, Mark Longair, Tobias Pietzsch, Stephan Preibisch, Curtis Rueden, Stephan Saalfeld, Benjamin Schmid, Jean-Yves Tinevez, Daniel James White, Volker Hartenstein, Kevin Eliceiri, Pavel Tomancak, and Albert Cardona. *Fiji: an open-source platform for biological-image analysis.* Nature Methods, 2012. <https://doi.org/10.1038%2Fnmeth.2019>
3. TrackMate v6.0.3 (Fiji plugin)
    * Jean-Yves Tinevez, Nick Perry, Johannes Schindelin, Genevieve M. Hoopes, Gregory D. Reynolds, Emmanuel Laplantine, Sebastian Y. Bednarek, Spencer L. Shorte, and Kevin W. Eliceiri. *TrackMate: An open and extensible platform for single-particle tracking*. Methods, 2017. <https://doi.org/10.1016/j.ymeth.2016.09.016>
    * I found the code examples for scripting TrackMate to be especially helpful. My tracking script is built from an amalgamation of these.
4. Matlab-ImageJ (Fiji update site)
    * <https://imagej.net/scripting/matlab>
5. CircStat
    * Berens, Philipp. *A MATLAB Toolbox for Circular Statistics.* Journal of Statistical Software (2009). DOI: [10.18637/jss.v031.i10](https://www.jstatsoft.org/article/view/v031i10)
6. Java Heap Cleaner
    * Davide Tabarelli. *Java Heap Cleaner.* MATLAB Central File Exchange (2013). URL: <https://www.mathworks.com/matlabcentral/fileexchange/36757-java-heap-cleaner>. Retrieved August 28, 2021.