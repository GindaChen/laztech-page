## VLASS data products

K.H. Yuen 						Sep 28 2018

#### NRAO archine

- [Old site]:(archive.nrao.edu)

- [New site]: (archive-new.nrao.edu): esp /vlass

- [Quick-look]: (https://archive-new.nrao.edu/vlass/quicklook/VLASS1.1/)

#### Data processing

- Sky tiers from $RA =-40\deg \sim90 \deg$ 

- Quality check:
  ![1538163652942](C:\Users\KHYuen\AppData\Roaming\Typora\typora-user-images\1538163652942.png)

- Quicklook: $1 \deg^2$ tiles 

  > Tile name,  image center J2000, 

- [Instruction]: (go.nrao.edu/vla-pipe) , [Pipeline guide]:(go.nrao.edu/vla-casa-tut)

- [Pipeline Weblog]: (casaguides.nrao.edu)

#### Data Access

- https://archive-new.nrao.edu/vlass/quicklook/ will have everything you need

- FIlename in this structure

  > - VLASS epoch (1.1=first half of first epoch; 1.2=second half of first epoch; 2.1=first half of second epoch, etc.)
  > - Product type (ql=Quick Look imaging)
  > - Tile name
  > - Image phase center
  > - Pixel size x10 (10=1.0arcsec)
  > - Bandwidth in MHz
  > - Version number
  > - Stokes type

- What we will look![1538165764857](C:\Users\KHYuen\AppData\Roaming\Typora\typora-user-images\1538165764857.png)

  > - `rms.subim.fits` : deviation image
  >
  > - `subim.fits ` :the actual image
  >
  > - `CASA_command.log` :pipeline command log
  >
  > - `casa_pipescript.py` :high-level CASA pipeline script
  >
  >   ``````python
  >   __rethrow_casa_exceptions = True
  >   context = h_init()
  >   context.set_state('ProjectSummary', 'proposal_code', 'VLA Prop Code')
  >   context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
  >   context.set_state('ProjectSummary', 'telescope', 'EVLA')
  >   context.set_state('ProjectSummary', 'piname', 'unknown')
  >   context.set_state('ProjectSummary', 'proposal_title', 'unknown')
  >   try:
  >       hifv_importdata(nocopy=True, vis=['VLASS1.1.sb34917571.eb34994250.58152.889074189814.ms'], session=['session_1'])
  >       
  >   	# KH: The brelow line reads the parameter list to change image size
  >       hif_editimlist(parameter_file='parameter.list') 
  >       # --
  >       
  >       hif_transformimagedata(datacolumn='corrected', modify_weights=False, clear_pointing=True)
  >       hif_makeimages(hm_masking='none', hm_cleaning='manual')
  >       hifv_pbcor(pipelinemode="automatic")
  >       hif_makermsimages(pipelinemode="automatic")
  >       hif_makecutoutimages(pipelinemode="automatic")
  >   finally:
  >       h_save()
  >   ``````
  >
  > - `parameter.?????.list` :Input parameter
  >
  >   ``````
  >   # Default values
  >   # imsize = [7290,7290]
  >   # serach_radius_arcsec = 1000
  >   # cycleniter = 500
  >   
  >   editmode='add'
  >   imagename='VLASS1.1.ql.T01t02.J003228-363000.10.2048.v1'
  >   phasecenter='J2000 00:32:28.328 -36.30.0.0000'
  >   ``````
  >
  > - `pipeline...xml`: CASA version used
  >
  > - `weblog.tgz`: weblog, can run offline
  >
  >   ![1538166269226](C:\Users\KHYuen\AppData\Roaming\Typora\typora-user-images\1538166271629.png)

- Other parameters you can modify:

  > Check weblog for keywords like: `hif_editimlist`, `hifv_editimlist`

- Want to run pipeline? Typical dates to work: 14 days.

  > Suggested way to do
  >
  > 1. Find your interested tiles in the weblog archive
  > 2. They will contain the pipeline 
  > 3. Everything is there with the numbers on the log hinting you what is going on

- Notes for downloading the data: 

  > - For more data $\rightarrow$ Go old NRAO archive.
  >
  > - The new NRAO page (archive-new) allows text search + filter.
  >
  > - VLASS **1**.2 : 2nd half of the sky, **1**st epoch
  > - Calibrated files can also be download from the sit
  > - Saver to download raw data, and convert the data using the measurement sets in `CASA`
  > - SDM = Science Data Mod