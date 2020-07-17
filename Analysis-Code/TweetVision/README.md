# TweetVision: Basic Use Instructions
Unlike TweetVisionLite, TweetVision has functions for electrophysiology data. It should be used for .DAT files. Tweet Vision has two major functions: it can be used to annotate audio recordings just like TweetVisionLite, but it can also be used to visiually check spike time files, manually edit them, and save them into cell files, which contain records of when different neurons are firing according to the spikes seen in the electrophysiology data. 
## Interface Description:
**UPDATE UPDATE UPDATE**
![TweetVision GUI](https://user-images.githubusercontent.com/18174572/41117110-1b2fc5f0-6a5a-11e8-8ebf-f75af9ef04a7.png)
* **DAT File Control**: The DAT file control panel, located in the top right corner of the application, is used to load folders and navigate the files they contain.
    * Load File: opens a finder window that allows you to choose the folder you want to annotate
    * Previous/Next: scroll through .DAT files in the loaded folder
    * List View: shows the names of all the files in the loaded folder
    * Quarentine/Empty Quarentine: a now outdated tool (due to the existence of SongBlaster), these buttons were used to flag/unflag certain files for deletion
    * Delete: delete quarentined files
* **Annotation**: The tools in this panel are used to create annotation files and syllable patterns (similar to those created in TweetVisionLite). 
    * New Annotation: creates a new annotation file; you must creat a new annotation for each folder you load before doing anything else; name according to conventions described in the use instructions
    * Load Annotation: opens a finder window where you can select a previously created annotation file to continue to edit
    * Save Template: Save the current syllable pattern
    * Add Syllable: press this button prior to clicking and dragging to select sections on the spectrogram to create syllables
    * Delete Syllable: click a selection on the spectrogram and then press this button in order to delete an unwanted syllable
    * View Template: Load/view a previously created syllable pattern
    * Show Audio Power: ?????????
    * Edge Finder: OBSOLETE
    * Delete Pause: joins two separate syllables
    * Load Annotation: load a previously generated annotation file to add to it
    * View Template: view a previously generated syllable pattern
    * Infused Drugs/Directed-Undirected: provide more information about the annotations you are creating
* **Audio and Spectrogram Displays**: These are identical to those in TweetVisionLite. Here is where you will segment and label syllables in the audio files. 
    * Audio Envelope: a graphical representation of the amplitude of the songs; click and drag in the audio envelope to zoom in on the spectrogram below
    * Spectrogram: a visual representation of the frequency and amplitude of the bird's calls; the x-axis represents time, the y-axis represents frequency, and the color denotes amplitude; add new syllables by clicking and dragging on this view and numberng your selections
* **Segmenting**: These controls assist in segmenting songs and assigning numbers to different syllables.
     * Multi-Segment: auto-segments entire file
    * Frame Segment: auto-segments all syllables in the currently selected box in the audio envelope
    * Clear All: removes segmentation on spectrogram
    * Reset Params: resets all parameters to their default values
    * Power Env/Dynamic/Seg Guide: segmentation settings
    * Syllables/Bouts: different parameters that can be adjusted to improve the accuracy of the auto-segmenting functions
* **E-Phys Data Display**: These controls are to adjust the electrophysiology graphs on the immediate left of the panel. You can choose to display the voltage, power, spectrogram, or spectrum of the current .DAT file.
    * Voltage/Power/Spectrogram/Spectrum: select which display you would like to view
    * CMS/Differential: alter E-Phys data displays
    * Scale Freq: scales the frequency of the spectrum graphs based on the number you enter in the box below
    * Pre-processing/Update: pre-processing opens a new window in which you can alter the filtering on the electrophysiology data, update applies these new filters to the data
    * Min/Max: sets the bounds of the y-axis in the electrophysiology data views
    * Sync Scales: ensures the scales are aligned
* **Spike Sorting Control**: This panel is used to manage different cell files that you either create or load in.
    * New Cell: create a new cell file
    * Load Cell: load in a previously created cell file
    * Unload File: unload a currently loaded cell file
    * Display Box: displays all currently loaded cell files
    * Show Spike Rasters: display spike markers on the axis below the e-phys displays
* **E-Phys Graphics & Controls**: These graphs and buttons are used to create and edit cell files, which allow you to cluster different spikes in the electrophysiology data as specific action potentials. Each of the four graphs corresond to a different channel of e-phys data. These controls are typically used to edit cell files that were generated by the Wave_clus package and verify their accuracy.  
    * CellFile dropdown: select which of the cell files from the spike sorting control panel to display
    * CellNum dropdown: select which cell number from the currently selected cell file to display
    * "+": click on axis and add a spike
    * "-": get rid of a spike
    * Broom: clear all spikes
    * Recycle: reload spikes from file, discard manual changes
    * "!": insert a spike whether or not a marker exists
* **Automated Syllable ID**: This panel (in the bottom right) allows you to automate the process of identifying the syllables in each audio file. In order to implement this functionality, some syllables must be annotated manually and a neural network must be trained.
    * Train Network: train a network to identify syllables for you
    * Load Network: load a previously trained network for use in this folder; calls should be fairly consistent for the same bird, so once you train one neural network for a particular bird you should be able to reuse it for any recordings from that bird
    * Save network: allows you to save a trained network for later use (adhere to naming conventions described in use instructions)
    * ID Syllables: use the loaded neural network to classify all segmented syllables in the current file
    * Batch Process: automates segmentation AND syllable ID's across multiple flies
    * Overhang/SylDurWeighting: neural network settings
    * Decision Pnt: confidence threshold for neural network assigning syllables or "unknown"
    * Time: assigns time frame for batch process
* **Annotation Progress**: this panel (located in the center on the bottom of the GUI) allows you to see how much you've annotated and how many of each pattern of syllables you have identified 
* **Main Control**: These are controls for saving and plotting data; they are not part of the main functionality of the application. 
    * Save Current Data: saves name, audio, neural channels, and filtering parameters of current file
    * Plot Correlation: plots the correlation between two channels of your choice

##  Instructions:
### Syllable Annotation
#### Manual Use
1. Load folder using file control panel, select a folder containing .wav files, and click open
2. An "Input Experiment Parameters" window will appear; all you have to do is enter the bird's name
3. A warning will appear saying that you have not yet created an annotation file; this is fine
4. Click "new annotation", name according to the convention Birdname_YYMMDD_annotation.mat
5. Double click a file name in the file control panel. This will bring up the audio envelope and spectrogram for that file.
6. Click and drag on the audio envelope to zoom in
7. Click **add syllable**, click and drag around syllable in spectrogram
8. Click on newly created syllable and then number it
9. Add syllables to the syllable pattern by doubleclicking them
10. To save a syllable pattern, use "save template" button in the annotation panel
11. Get annotating! :)

* Note that annotated files turn green after you've worked with them to indicate that they have already been annotated! 

#### Automated Use
First, follow steps 1-5 in the manual use instructions. 

6. Click and drag on the audio envelope to select a subset of the recording. Then click "frame segment" to auto-segment this subset into syllables. Alternatively, use "multi-segment" to auto-segment the entire file at once
7. If you have done enough manual IDing of different syllables (preferably 100 examples of each syllable), you can train a neural network to do the rest for you! Click the "Train Network" button in the Automated Syllable ID panel
8. Once the "Train Network" button turns green, click the "ID syllables" button to automatically classify each syllable in the selected file
----------------------------
The **Batch Process** can be used to automatically segment and annotate multiple files at once. Note that a network must be trained for this function to be used. Select the start and end times for the files that you want annotated and then press the button. 

### Electrophysiology and Cell Files
View and manually edit the spiketime files. 
1. Create a new cell using the "New Cell" button in the spike sorting control panel or load a previously created cell with the "load cell" button and check the "Show Spike Rasters" box
2. View your desired cell file using the CellFile dropdown menu next to the graph of a channel of electrophysiology data
3. Now you can edit the spikes however you would like: click anywhere in the lower display and then use the "+", "X", broom, recycling, and "!" buttons as described above until the spikes are to your liking
4. When you are satisfied with the spike time file, hit "next" in the DAT File Control panel to save the changes to the current cell file
