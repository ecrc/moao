#!/bin/sh

datafile="../datafile/check/"


maxrefine=1
maxobs=1
snap_per_night=8

./test_moao_lapack --filePath="$datafile/nx456_nLayers2_wfs6_Nssp10/"    --suffix="_nx456_nLayers2_wfs6_Nssp10"    --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night > log_val 2> err_val
./test_moao_lapack --filePath="$datafile/nx456_nLayers10_wfs6_Nssp10/"   --suffix="_nx456_nLayers10_wfs6_Nssp10"   --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night>>log_val 2>>err_val 
./test_moao_lapack --filePath="$datafile/nx5120_nLayers2_wfs8_Nssp29/"   --suffix="_nx5120_nLayers2_wfs8_Nssp29"   --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night>>log_val 2>>err_val 
./test_moao_lapack --filePath="$datafile/nx5120_nLayers10_wfs8_Nssp29/"  --suffix="_nx5120_nLayers10_wfs8_Nssp29"  --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night>>log_val 2>>err_val 
./test_moao_lapack --filePath="$datafile/nx10112_nLayers2_wfs8_Nssp41/"  --suffix="_nx10112_nLayers2_wfs8_Nssp41"  --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night>>log_val 2>>err_val 
./test_moao_lapack --filePath="$datafile/nx10112_nLayers10_wfs8_Nssp41/" --suffix="_nx10112_nLayers10_wfs8_Nssp41" --maxrefine=$maxrefine --maxobs=$maxobs  --snap_per_night=$snap_per_night>>log_val 2>>err_val 
