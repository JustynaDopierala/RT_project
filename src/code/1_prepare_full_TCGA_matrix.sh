#Author: Justyna Dopierala
#Date: Nov 2023

###########################################################################
### Combine TPM columns from multiple files
###########################################################################

# read TCGA data input path
dir=$1
counter=0
for f in "$dir"/*;do
  counter=$[counter +1]     # count files
  a="$(basename -- $f)"     # extract file name
  header="${a}_TPM"         # make new header
  if [ $counter = 1 ]       # if first file
  then
    awk '{print $1,$2,$3}' $a > collapsed_output.tsv      # print the first 3 columns to another file (annotations)
    awk '{print $7}' $a > temp.tsv                        # print the column with TPM data to a temp file
    sed -i "1s/.*/$header/" temp.tsv                      # replace old header for new header
    paste -d'\t' collapsed_output.tsv temp.tsv > tmpout && mv tmpout collapsed_output.tsv     # join annotations with TPM
    rm temp.tsv                                                                               # remove temp file
  else                                # if not first file
    awk '{print $7}' $a > temp.tsv    # print the column with TPM data to a temp file
    sed -i "1s/.*/$header/" temp.tsv  # replace old header for new header
    paste -d'\t' collapsed_output.tsv temp.tsv > tmpout && mv tmpout collapsed_output.tsv     # join annotations with TPM
    rm temp.tsv                       # remove temp file
  fi
done

output_dir=$2
cp collapsed_output.tsv $output_dir   # copy to a destination directory
