# SongBlaster: Basic Use Instructions
Songblaster is used to scrub audio recording data of background noise. It receives a folder of audio recordings from one bird during one specific day.

Note: These instructions describe the *typical* use of the SongBlaster package.

## Interface Description
![songblastergui](https://user-images.githubusercontent.com/18174572/40849024-dc52fbe2-658e-11e8-8578-4a9f1ff0a113.png)

* **File Control**: The file control panel (located in the upper right hand corner of the SongBlaster GUI) allows you to load folders containing audio files and delete those deemed to be background noise. It has several features:
    * Top bar: the top black bar displays the name of the folder that you are working on
    * Load Folder: this button opens up a finder window that allows you to chose the folder you want to scrub
        * Load Folder Checkboxes: allow you to select the specific data types you are interested in
    * Auto Trash: Allows you to select all audio recordings that are shorter than a specified recording (entered in the box above)
    * Reset Trash: Unchecks all of the files currently marked for deletion.
    * Show Trash Bin: sets audio display to only show files that have been marked for deletion
    * Delete: deletes all currently selected files (!!Notice: these files cannot be recovered!!)
* **Navigation**: Located in the box to the left of File Control, the navigation tools allow you to scroll through the recordings and change what is shown in the waveform display. The left scroll bar changes the page number and the right scroll bar changes the particular file that is displayed. The application will display the waveform and spectrogram of whichever particular audiofile is selected in the navigation box in the bottom right corner of the GUI.
* **Stay/Go**: In this box you select all of the files that you want to delete. As you scroll to different pages using the navigation bar, it remembers files you have selected on other pages.
* **Waveform Display**: The entire left side of the SongBlaster GUI is dedicated to displaying the waveforms of the audio/electrophysiology files in the loaded folder. It shows 10 files at a time. The bottom right corner gives a more detailed display of one particular file (shows a spectrogram and a zoomed-in waveform). 
* **Data Display**: 
    * Left Scroll Bar: this bar allows you to adjust the amplitudes of the waveforms so they are more meaningful; the scroll bar in the bottom right next to the zoomed-in waveform does the same thing for that particular graph
    * Dropdown Menu: not important for .WAV files, but for .DAT files this tool allows you to choose the channel shown in the waveform display; it can show audio or any of the electrophysiology channels contained in the file; altering this dropdown menu also affects the spectrogram/large waveform display in the bottom right
    * CMS Checkbox: creates an average reading from the electrophysiology data channels and subtracts it from all of them; results are shown in the bottom right detailed view
    * Vref Checkbox and Dropdown Menu: Allows you to select one of the electrophysiology data channels (via the dropdown menu) and subtract it from the rest; results are shown in the bottom right detailed view
    * Preprocessing/Update: Not important
* **Distribution Histogram**: The distribution histogram gives the user a picture of when during the day the different audio files were recorded. The x axis denotes time of day and the y axis shows the number of recordings taken at that specific time.

## Instructions: 
1. Click the Load Folder button in the File Control panel. This will open a finder window. Navigate to the folder containing your data files and click open.
2. Using the navigation panel, scroll through the audio files to determine which ones are background noise. Check the garbage files off as you go. If you prefer, use the auto trash feature by entering the number of a recording (e.g. for the third entry on page two enter 23) into the box above the "Auto Trash" button and clicking it. This will select all of the recordings shorter than the one you selected. 
3. Once you have selected all of the files that you want to get rid of, press the "delete" button. This will **PERMANENTLY** delete all of the files you have selected, so be careful!
4. The rest of the functions of this application are just used to view the files (to judge if they are background noise or actual birdsong). Use as necessary. 
