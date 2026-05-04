#!/bin/bash

# List of all mass points
mass_points=(15 20 25 30 35 40 45 50 55 60 62)

# Loop and call ROOT macro with each mass point
for mA in "${mass_points[@]}"; do
  echo "Running for mA = $mA"
  root -l -b -q "H2AA2bbgg_analysis_treeread_18_for_pas_t2.C($mA)"
  echo "Done for mA = $mA"
done

