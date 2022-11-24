## MFT Apply Survey Align Parameters 

This repo contains the raw ideal/survey positions used to compute alignment parameters for the MFT (at two levels: chips, disks)

It contains:
    RunSurveyAlignWorkflow.cxx
    ReadMFTSurvey.C
    ModuleTransform.h
    ModuleTransform.cxx
    rootlogon.C
    
geometry files needed:
    Geometry.xml
    o2sim_geometry.root
    
and mftalign folder with positionsHIC (all the survey for chips)
other folders: 
    temp_files
    data_txt_files
    

Get the survey with: 
    root -l -b ReadMFTSurvey.C
    
To be run with: 
    root -l -b RunSurveyAlignWorkflow.cxx
