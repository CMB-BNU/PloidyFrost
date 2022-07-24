thread=1
dir="."
name="sample"
readsidlist="idlist"
suffix=".fastq.gz"
script=""

time sh $script/1.trim $thread $dir $readsidlist $suffix 
time sh $script/2.kmc_db $thread $dir $name $readsidlist $suffix
coverageL=$(PloidyFrost cutoffL ${dir}/hist_${name})
coverageU=$(PloidyFrost cutoffU ${dir}/hist_${name})
time sh $script/3.filter $thread $dir $name $coverageL
time sh $script/4.bifrost $thread $dir $name
time sh $script/5.ploidy $thread $dir $name $coverageL $coverageU

