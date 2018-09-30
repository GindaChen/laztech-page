# VLA Skysurvey Talk

`@Talker: Amy, 9/28 2:30-3:45pm`



## Summary

### Useful Links

#### Quicklook & Downloads

- **[Scrollable Selection Page](http://archive-new.nrao.edu/vlass/HiPS/VLASS1.1/Quicklook/)**: You can actually play with the glob, zoom in to see high resolution image, and find where to download.
- **[Quick look](https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/)**
- **[VLASS Resource New Archive Page](https://archive-new.nrao.edu/vlass)**: The main page including all the information about observation: [Summary of the cubes](https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php), [Calibration and pipeline](https://archive-new.nrao.edu/vlass/weblog/), [Product - Quicklook images](https://archive-new.nrao.edu/vlass/quicklook/), and of course [the scrollable selection](http://archive-new.nrao.edu/vlass/HiPS/).
- [**Instruction for installation and operation**](https://archive-new.nrao.edu/vlass/weblog/quicklook/): 

#### General Information

- Old Archive: https://archive.nrao.edu
- New Archive: https://archive-new.nrao.edu/

#### Observation Details 

- [Pipline guide](https://archive-new.nrao.edu/vlass/weblog/quicklook/)



### Quickstart 

**Our goal: to get a calibrated dataset.** The most important thing for us is to get a **calibrated dataset** at a specific zone. Hence we should find the **single epoch product**, which has not been released yet.



Step 1: Find your cube in the Quicklook page or [Summary page](https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php) where the basics of data is described as below

```csv
Tile     Dec min  Dec max   RA min   RA max    Observing     Observation         Quick Look
Name      (deg)    (deg)     (h)      (h)      Epoch         Date (UT)           Status
-------------------------------------------------------------------------------------------------
T01t01   -40.00   -36.00     0.00     0.50     VLASS1.1      2018-02-07               100% imaged
T01t02   -40.00   -36.00     0.50     1.00     VLASS1.1      2018-02-03               100% imaged
...
```



Step2: 



## Basics

### Data Types

There are several different datatypes:

- Raw data
- Calibrated
- ...
- **Single epoch product** : the final product

### Epochs

- v1.1: we have observed 
- v1.2: will begin by 2019/01

### Data Processing

#### 

- Sky tiers from $RA =-40\deg \sim90 \deg$ 

- Quality check:
  ![1538163652942](C:\Users\KHYuen\AppData\Roaming\Typora\typora-user-images\1538163652942.png)

- Quicklook: $1 \deg^2$ tiles 

  > Tile name,  image center J2000, 

- [Instruction]: (go.nrao.edu/vla-pipe) , [Pipeline guide]:(go.nrao.edu/vla-casa-tut)

- [Pipeline Weblog]: (casaguides.nrao.edu)



### File Structure

FIlename in this structure

> - VLASS epoch (1.1=first half of first epoch; 1.2=second half of first epoch; 2.1=first half of second epoch, etc.)
> - Product type (ql=Quick Look imaging)
> - Tile name
> - Image phase center
> - Pixel size x10 (10=1.0arcsec)
> - Bandwidth in MHz
> - Version number
> - Stokes type



## What data do we have

**Calibration and Imaging Weblogs**

Calibration data (deal with satelites, and other calibration dataset )

*Two types of calib*: Oppression and that not-oppresion 



[40 per tile](https://archive-new.nrao.edu/vlass/weblog/quicklook/)

```
Tilename, Image center (J2000 RA/sec), Unique Pipline Job id
T, J, P
```





**Pipeline Overview** (I have no idea what it is...)

standard observing

(? too detail... underflag? overflag?)





**VLA casa pipeline (detailed)**



**Idea:** 

1. Pipe line is created real-time
2. Diagonostic is updated afterwards

Imprtant steps:

- Apply calibration
- Flag data 
- Generate new dataset
- Diagnostic, and QA





### Example



#### Overview of the pages functionality



https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t2-1.html?sidebar=sidebar_TSKY0001_sb32429825_eb32434215_57582_960989803236_ms&subpage=t2-1_details.html

- [Overview of the pipeline](https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t1-1.html)
- [By Task](https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t1-4.html)



#### Diagonostic: Plot summary (item #19)

https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t2-4m.html?sidebar=sidebar_stage19&ms=all&subpage=t2-4m_details.html



**Point source and not point source**: if flag and sharp, it's point source

**Polarized & Unpolarized**: straight line vs scatter



#### Details

[3. VLA Deterministic FlaggingBACK](https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t2-4m.html?sidebar=sidebar_stage3&ms=all&subpage=t2-4m_details.html)



![image-20180928145817214](assets/image-20180928145817214.png)





![image-20180928145743974](assets/image-20180928145743974.png)

**More than 20% data flaged => data might be wrong**



[5. Prior calibrationsBACK](https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t2-4m.html?sidebar=sidebar_stage5&ms=all&subpage=t2-4m_details.html)



https://archive-new.nrao.edu/vlass/weblog/calibration/VLASS1.1_GOODSN_P27v1_2018_03_15T22_51_17.728/pipeline-20180316T000655/html/t2-4m.html?sidebar=sidebar_stage5&ms=all&subpage=spgain-TSKY0001.sb32429825.eb32434215.57582.960989803236.ms.html

Unstable data: peak/valley

Compress data: drop

![image-20180928150022404](assets/image-20180928150022404.png)



















## Main theme: Quickimaging



Final images: [the cutout images](https://archive-new.nrao.edu/vlass/weblog/quicklook/VLASS1.1_T01t01.J000228-363000_P25905v1_2018_04_12T23_15_12.148/pipeline-20180412T231946/html/t2-4m.html?sidebar=sidebar_stage7&ms=all&subpage=t2-4m_details.html)



- subimage (beam energy distribution)
- residual image:







## VLA Data Access: finally

[Quicklook: where you grab images](https://archive-new.nrao.edu/vlass/quicklook/)

![image-20180928151339459](assets/image-20180928151339459.png)

- VLASS epoch (1.1=first half of first epoch; 1.2=second half of first epoch; 2.1=first half of second epoch, etc.)
- Product type (ql=Quick Look imaging)
- Tile name
- Image phase center
- Pixel size x10 (10=1.0arcsec)
- Bandwidth in MHz 
- Version number
- Stokes type



[All the images of VLASS1.1](https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/)

- `rms.subim.fits` : deviation image
- `subim.fits ` : the actual image (contain noise, $1\times1\  deg^2$ image)
- `CASA_command.log` : pipeline command log
- `casa_pipescript.py` : high-level CASA pipeline script
- `parameter.list` : Parameter list for your own pipeline generation
- `pipeline_manifest.xml`: CASA version used
- `weblog.tgz`: weblog, can run offline



 [**Example Page**](https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/T01t01/VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1/)

- `VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1.I.iter1.image.pbcor.tt0.rms.subim.fits 16-Apr-2018 15:43   53M `
- `VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits `: 

```
3. casa_commands.log
4. casa_pipescript.py
5. parameter.list
6. pipeline_manifest.xml
7. weblog.tgz
```



1. deviation image
2. actual image (1x1 deg^2 image)
3. & 4.-> we can run the pipeline ourself.

### Generate your own pipeline

[casa_pipescript.py](https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/T01t01/VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1/casa_pipescript.py)

```python
__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'proposal_code', 'VLA Prop Code')
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
context.set_state('ProjectSummary', 'piname', 'unknown')
context.set_state('ProjectSummary', 'proposal_title', 'unknown')
try:
    hifv_importdata(nocopy=True, vis=['VLASS1.1.sb34916486.eb35006898.58156.86241219907.ms'], session=['session_1'])
    hif_editimlist(parameter_file='parameter.list')
    hif_transformimagedata(datacolumn='corrected', modify_weights=False, clear_pointing=True)
    hif_makeimages(hm_masking='none', hm_cleaning='manual')
    hifv_pbcor(pipelinemode="automatic")
    hif_makermsimages(pipelinemode="automatic")
    hif_makecutoutimages(pipelinemode="automatic")
finally:
    h_save()
```



[Paramlist](https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/T01t01/VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1/parameter.list)

```python
editmode='add'
imagename='VLASS1.1.ql.T01t01.J000228-363000.10.2048.v1'
phasecenter='J2000 00:02:28.328 -36.30.0.0000'

####### Change these parameters to generate more acurate 
search_radius_arcsec= 1000
imgsize=[7290, 7290]
cycleniter= 500 # cycle iteraton
```





(??? How to calibrate the search_radius to the search area?)



## How to change the param so that we have smaller images

![Image from iOS](assets/Image from iOS.jpg)







## Raw data: 

**This is really big, really raw, not calibrated, and it's not good to run alone :)**

[Calibration/compression/a_specific_pipeline/html/the_ms_file.ms](Calibration/compression/a_specific_pipeline/html/the_ms_file.ms)





## New Archive support only a part of data; more are archived in the old page

archive.nrao.edu

archive-new.nrao.edu

| -- [Project view](https://archive-new.nrao.edu/archiveIface/#/?pageNumber=1&pageSize=25&sortOnColumn=obs_stop&sortDir=desc&currentResultsView=project&showFilters=false)

| -- [Observation view](https://archive-new.nrao.edu/archiveIface/#/?pageNumber=1&pageSize=25&sortOnColumn=obs_stop&sortDir=desc&currentResultsView=classic&showFilters=false)







execution block id -> execution block



# Downloads:

![image-20180928153023521](assets/image-20180928153023521.png)

![image-20180928153032157](assets/image-20180928153032157.png)



download only calibrator scan ...





[Choose download data format:](https://archive-new.nrao.edu/archiveIface/)

- SDM tables only (no visibilities) = metadata only
- SDM-BDF dataset (all files) = All tables + visibilities, **but have to run a calibration**
- Basic Measurement Set (uncalibrated) = CASA data format
- **Calibrated Measurement Set = Not fully tested, but the most calibrated data set**

Look at weblog first, see if the clibration good, and then play with it, and then download calibrated ms.

(ms = measurement set)





You can choose the casa that's going to run

![image-20180928153327888](assets/image-20180928153327888.png)







Should we apply telescope flags -- recommended!

average in time? (beaware when averaging uncalibrated ms)





![image-20180928153432358](assets/image-20180928153432358.png)









[Events](https://go.nrao.edu/vla-we)





# Next Cycle: Start at 2019/01
