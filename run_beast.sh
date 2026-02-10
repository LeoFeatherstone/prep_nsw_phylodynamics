arr=(re_averaged.xml cluster_161.xml cluster_168.xml cluster_37.xml dpp.xml)
for xml in "${arr[@]}"; do
    beast2 -overwrite -seed 4321 -beagle xml/$xml
done