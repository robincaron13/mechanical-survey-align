## MFT Apply Survey Align Parameters 

This repo contains the raw ideal/survey positions used to compute alignment parameters for the MFT (at two levels: chips, disks)

It contains:
    `RunSurveyAlignWorkflow.cxx`,
    `ReadMFTSurvey.C`,
    `ModuleTransform.h`,
    `ModuleTransform.cxx`,
    `rootlogon.C`
    
geometry files needed:
    `Geometry.xml`,
    `o2sim_geometry.root`
    
and `mftalign/positionsHIC` folder with (all the survey txt data for chips on HICs)
other folders: 
    `temp_files`,
    `data_txt_files`
    
Enter in O2 with: `alienv enter O2/latest-master-o2`

Get the survey with: 
    `root -l ReadMFTSurvey.C`
    
To be run with: 
    `root -l RunSurveyAlignWorkflow.cxx`
