ckfPropSteps="100000000"
#kfPropSteps="10000"
maxSeeds="10000"
ckfchi2max=30
#default is 2.74
cotThetaMax=3.0
#default is 0.1 
radLengthPerSeed=0.03

matInputFile=/mnt/c/Users/lhaza/Desktop/stcf-tgeo/mat-STCF_tracker_alt5_manual_MDI.json
geoTGeoJsonFile=/mnt/c/Users/lhaza/Desktop/stcf-tgeo/tgeo_STCF_tracker_config_alt.json
geoTGeoRootFile=/mnt/c/Users/lhaza/Desktop/stcf-tgeo/STCF_tracker.root

inputDir=/mnt/c/Users/lhaza/Desktop/run_data
#inputFile1=fullsim_2e2pi.root 
#inputFile2=fullrec_2e2pi.root
inputFile1=muon_simulation.root
inputFile2=muon_rec.root
#inputFile2=fullrec_lambdalambdabar.root


outputDir=/home/lihao666/acts-core/build/bin/run/$tag

#--ckf-initial-variance-inflation=1:1:1000:1000:1000:1
nbranch=1
#--input-files=${inputFile}
run="/home/lihao666/acts-core/build/bin/ActsExampleCKFTracksSTCF_houghseed 
--seed-cottheta-max $cotThetaMax 
--seed-rad-length-per-seed $radLengthPerSeed 
--seed-max-seeds $maxSeeds 
--ckf-prop-steps $ckfPropSteps  
--ckf-selection-nmax $nbranch 
--ckf-selection-chi2max=${ckfchi2max} 
--ckf-initial-variance-inflation=1:1:1000:1000:1000:1 
-j 1 
--mat-input-type file 
--mat-input-file=${matInputFile} 
--geo-tgeo-filename=${geoTGeoRootFile}  
--geo-tgeo-jsonconfig=${geoTGeoJsonFile}  
--input-dir=${inputDir} 
--input-files=${inputFile1}
--input-files=${inputFile2}
--output-dir=${outputDir} 
--bf-constant-tesla=0:0:1  
--geo-selection-config-file=/mnt/c/Users/lhaza/Desktop/stcf-tgeo/geoSelection-STCFDetector.json"
eval $run
