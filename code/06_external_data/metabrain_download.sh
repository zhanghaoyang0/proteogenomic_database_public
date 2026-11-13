
# https://download.metabrain.nl/files.html

wget -c https://download.metabrain.nl/2021-07-23-release/columnDescription.txt

for trait in basalganglia cerebellum cortex hippocampus spinalcord; do
for chr in {1..22}; do
if [ "$trait" == "cerebellum" ]; then
  npcs=60
elif [ "$trait" == "cortex" ]; then
  npcs=80
elif [ "$trait" == "spinalcord" ]; then
  npcs=20
else
  npcs=30
fi
wget -c https://download.metabrain.nl/2021-07-23-release/2021-07-23-${trait}-EUR-${npcs}PCs-chr${chr}.txt.gz 
wget -c https://download.metabrain.nl/2021-07-23-release/2021-07-23-${trait}-EUR-${npcs}PCs-chr${chr}.txt.gz.tbi
done
wget -c https://download.metabrain.nl/2021-07-23-release/2021-07-23-${trait}-EUR-${npcs}PCs-TopEffects.txt.gz 
done