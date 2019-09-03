
cd archive/output000

rm *.nc

rm *profile*.*ta

archive_path='/g/data/v45/jm5970/particle_release/30d'
experiment_path={0}

mkdir '$archive_path/$experiment_path'

rsync -qar --remove-source-files ./* '$archive_path/$experiment_path'


