

cd sub-28680644/091008_00_28680644
datalad remove *
rm -r *
cp -r /home/gabridele/backup/ROSMAP/batch_2/BNK/090211/091008_00_28680644/3D* .

cp -r /home/gabridele/backup/ROSMAP/batch_2/BNK/090211/091008_00_28680644/Obl* . 

cp -r /home/gabridele/backup/ROSMAP/batch_2/BNK/090211/091008_00_28680644/P16* .

cp -r /home/gabridele/backup/ROSMAP/batch_2/BNK/090211/091008_00_28680644/rfMRI* .
cd ..
datalad save -m 'recopied sourcedata into folder, unexpected missing files'