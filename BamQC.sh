

if [ "$(basename $SHELL)" = "zsh" ]
then
  source ~/.zshrc
else
  source ~/.bashrc
fi

# $1 is the bam file, $2 is BED 12 file

# Check RSeQC installation
if [ -x "$(command -v bam_stat.py)" ]
then
  rm -f $(dirname $1)/bam_stat.txt
  bam_stat.py -i $1 > $(dirname $1)/bam_stat.txt
  infer_experiment.py -r $2 -i $1 >> $(dirname $1)/bam_stat.txt
else
  conda activate rseqc
  if [ $? -eq 1 ]
  then
    echo "Please install RSeQC by pip!"
  else
    rm -f $(dirname $1)/bam_stat.txt
    bam_stat.py -i $1 > $(dirname $1)/bam_stat.txt
    infer_experiment.py -r $2 -i $1 >> $(dirname $1)/bam_stat.txt
  fi
fi
