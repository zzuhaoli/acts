#!/bin/bash

#tag=CKF.estimated.chi2Cut15.maxPropSteps330.maxSeeds2.realmass
#tag=CKF.estimated.chi2Cut15.maxPropSteps10000.maxSeeds2.realmass
tag=CKF.estimated.PreditedDriftSign.chi2Cut30.maxPropSteps330.NoStepAjustError.maxSeeds2.realmass
#tag=CKF.estimated.chi2Cut15.maxPropSteps330.NoStepAjustError.LandauEloss.maxSeeds2.realmass
#tag=CKF.estimated.chi2Cut15.maxPropSteps330.NoStepAjustError.BetheEloss.maxSeeds2.realmass

momentums=("50" "75" "100" "125" "150" "175" "200" "250" "300" "350" "400" "450" "600" "800" "1000" "1200" "1400" "1600" "1800")
#momentums=("100")
#momentums=("50" "75")
absThetas=("36")
#absThetas=("90" "60" "30" "40")
particles=("mu-" "pi-")
#particles=("mu-")
#particles=("proton")

matInputFile=/home/xiaocong/Software/Oscar/acts/build/bin/mat-STCF_tracker_alt5_manual_MDI.json
geoTGeoJsonFile=/home/xiaocong/Software/Oscar/acts/build/bin/tgeo_STCF_tracker_config_alt.json
geoTGeoRootFile=/home/xiaocong/Software/Oscar/acts/build/bin/STCF_tracker.root
input=/nfs/xiaocong/Oscar_workspace/sim_momentum_angle_scan/pt/
ckfPropSteps="330"
maxSeeds="2"

for ((i=0; i<${#momentums[@]};++i)); do
for ((j=0; j<${#absThetas[@]};++j)); do
for ((k=0; k<${#particles[@]};++k)); do
momentum=${momentums[i]}
absTheta=${absThetas[j]}
particle=${particles[k]}
inputDir=${input}${particle}
inputFile=sim_${particle}_absThetaDeg_${absTheta}_${absTheta}_momentumPtMev_${momentum}_${momentum}.root
outputDir=/nfs/yi/acts_workspace/reco/scan/${tag}/${particle}/absThetaDeg_${absTheta}_${absTheta}_momentumPtMev_${momentum}_${momentum}

if [[ $particle == "proton" ]];then
   if [[ $momentum == "50" || $momentum == "75" || $momentum == "100" || $momentum == "125" || $momentum == "150" ]]; then 
     continue
   fi
fi

ckfchi2max=30
seedSigmaScattering="200"
ckfPropTolerance="0.0001"
ckfPropMass="139.57018"
if [[ $momentum == "50" || $momentum == "75" || $momentum == "100" ]];then
  seedSigmaScattering=2
fi

if [[ $momentum == "50" ]]; then
  ckfPropTolerance="0.01"
elif [[ $momentum == "75" || $momentum == "100" || $momentum == "125" ]];then
  ckfPropTolerance="0.0008"
fi

if [[ $particle == "mu-" ]];then
  ckfPropMass="105.6583755"
elif [[ $particle == "proton" ]];then
  ckfPropMass="938.272088"
fi


#Truth fitting
#KF.truthSmeared
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --smear-initial-variance-inflation=1000:1000:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1

#KF.estimated
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --fit-truth-estimated-seeds --smear-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json

#KF.estimated.chi2Cut15
#/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleTruthTracksSTCF --fit-chi2-cut=15 --fit-truth-estimated-seeds --smear-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json

#CKF

#CKF.estimated.chi2Cut15
run="/home/xiaocong/Software/Oscar/acts/build/bin/ActsExampleCKFTracksSTCF  --seed-max-seeds $maxSeeds --ckf-prop-steps $ckfPropSteps --ckf-prop-mass $ckfPropMass --ckf-selection-nmax 1 --ckf-selection-chi2max=${ckfchi2max} --ckf-initial-variance-inflation=1:1:1000:1000:1000:1 -j 1  --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1 --geo-selection-config-file=/home/xiaocong/Software/Oscar/acts/Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json"

echo $run
echo ======================================== 

eval $run

done
done
done




#2mu2pi
#inputFile=fullsim_2mu2pi_test.root
#outputDir=data/reco_STCF_OscarSim_reco_seeds/landauLoss/2mu2pi_test
#./ActsExampleCKFTracksSTCF --ckf-selection-nmax 1 --ckf-initial-variance-inflation=1:1:100:100:100:1 -j 1 --mat-input-type file --mat-input-file=${matInputFile} --geo-tgeo-filename=${geoTGeoRootFile}  --geo-tgeo-jsonconfig=${geoTGeoJsonFile}  --input-dir=${inputDir} --input-files=${inputFile} --output-dir=${outputDir} --bf-constant-tesla=0:0:1  --geo-selection-config-file=../../Examples/Algorithms/TrackFinding/share/geoSelection-STCFDetector.json -l 0 -n 50  > log 
