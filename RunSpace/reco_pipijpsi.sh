#!/bin/bash

#tag=v1.testLandau0p95.MDCDG
#tag=v1.1.testLandau0p95
#tag=v2.5seeds
#tag=v2.allseeds
#tag=CKF.estimated.chi2Cut15.maxPropSteps330.NoStepAjustError.betheEloss.5seeds
tag=CKF.estimated.PreditedDriftSign.chi2Cut50.maxPropSteps330.NoStepAjustError.allseeds


#Truth fitting
#v0
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --smear-initial-variance-inflation=1000:1000:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1

#v1
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --fit-truth-estimated-seeds --smear-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json

#v1.1
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --fit-chi2-cut=15 --fit-truth-estimated-seeds --smear-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json

#CKF

#run="/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleCKFTracksSTCF  --seed-max-seeds $maxSeeds --ckf-prop-steps $ckfPropSteps --ckf-selection-nmax 1 --ckf-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json"
#echo $run

#eval $run




#2mu2pi

ckfPropSteps="330"
maxSeeds="10000"
ckfchi2max=50

matInputFile=/home/xiaocong/Software/Oscar/acts/build/bin/mat-STCF_tracker_alt5_manual_MDI.json
geoTGeoJsonFile=/home/xiaocong/Software/Oscar/acts/build/bin/tgeo_STCF_tracker_config_alt.json
geoTGeoRootFile=/home/xiaocong/Software/Oscar/acts/build/bin/STCF_tracker.root

inputDir=/home/xiaocong/Software/Oscar/offline/Examples/JobSubmission/3RPC7PS/pipijpsi/
inputFile=fullsim.root
outputDir=/nfs/xiaocong/acts_workspace/reco/pipijpsi/$tag

run="/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleCKFTracksSTCF --seed-max-seeds $maxSeeds --ckf-prop-steps $ckfPropSteps  --ckf-selection-nmax 1 --ckf-selection-chi2max=${ckfchi2max} --ckf-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1 --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1  --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json"
echo $run
#-l 0 -n 50  > log 
